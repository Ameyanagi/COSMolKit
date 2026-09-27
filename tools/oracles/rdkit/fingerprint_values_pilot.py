"""Reference adapter only; Rust owns task selection, inputs and comparison."""
import json
import sys
import struct
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdDepictor

expected_version = sys.argv[1]
if rdBase.rdkitVersion != expected_version:
    raise RuntimeError(f"RDKit {rdBase.rdkitVersion} != {expected_version}")

def topology(mol):
    return {"atoms": [{
        "id": a.GetIdx(), "atomic_number": a.GetAtomicNum(),
        "isotope": a.GetIsotope(), "charge": a.GetFormalCharge(),
        "explicit_h": a.GetNumExplicitHs(), "no_implicit": a.GetNoImplicit(),
        "radicals": a.GetNumRadicalElectrons(), "aromatic": a.GetIsAromatic(),
        "hybridization": int(a.GetHybridization()), "chiral": int(a.GetChiralTag()),
        "permutation": a.GetUnsignedProp("_chiralPermutation") if a.HasProp("_chiralPermutation") else 0,
        "atom_map": a.GetAtomMapNum(),
    } for a in mol.GetAtoms()], "bonds": [{
        "id": b.GetIdx(), "begin": b.GetBeginAtomIdx(), "end": b.GetEndAtomIdx(),
        "order": int(b.GetBondType()), "aromatic": b.GetIsAromatic(),
        "conjugated": b.GetIsConjugated(), "direction": int(b.GetBondDir()),
        "stereo": int(b.GetStereo()), "stereo_atoms": list(b.GetStereoAtoms()),
    } for b in mol.GetBonds()]}


def molecular(row):
    profile = row["profile"]
    name, options = (profile, {}) if isinstance(profile, str) else next(iter(profile.items()))
    stage = "Parse"
    try:
        params = Chem.SmilesParserParams()
        params.sanitize = options["sanitize"] if name == "SmilesRead" else name != "SanitizeAll"
        params.removeHs = options["remove_hydrogens"] if name == "SmilesRead" else name != "SanitizeAll"
        params.allowCXSMILES = True
        params.strictCXSMILES = True
        params.parseName = True
        # The pinned Python v1 parser does not expose v2 skipCleanup;
        # its C++ conversion leaves that field at the default false.
        mol = Chem.MolFromSmiles(row["case"]["smiles"], params)
        if mol is None:
            raise ValueError("RDKit MolFromSmiles returned None")
        stage = "Operation"
        if name == "DistanceMatrix":
            matrix = Chem.GetDistanceMatrix(mol, useBO=options["use_bond_order"],
                useAtomWts=options["use_atom_weights"], force=True)
            return {"Matrix": {"dimension": mol.GetNumAtoms(), "values_bits": [
                struct.unpack(">Q", struct.pack(">d", float(value)))[0]
                for value in matrix.flat]}}
        if name == "MolecularWeight":
            value = Descriptors.MolWt(mol, onlyHeavy=options["only_heavy"])
            return {"Float64Bits": struct.unpack(">Q", struct.pack(">d", value))[0]}
        if name == "ExactMolecularWeight":
            value = Descriptors.ExactMolWt(mol, onlyHeavy=options["only_heavy"])
            return {"Float64Bits": struct.unpack(">Q", struct.pack(">d", value))[0]}
        if name == "MolecularFormula":
            return {"Text": rdMolDescriptors.CalcMolFormula(mol,
                separateIsotopes=options["separate_isotopes"],
                abbreviateHIsotopes=options["abbreviate_h_isotopes"])}
        if name == "SanitizeAll":
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        elif name == "Kekulize":
            Chem.Kekulize(mol, clearAromaticFlags=options["clear_aromatic_flags"])
        elif name == "AddHydrogens":
            params = Chem.AddHsParameters()
            params.explicitOnly = options["explicit_only"]
            params.addCoords = False
            params.addResidueInfo = False
            params.skipQueries = False
            mol = Chem.AddHs(mol, params)
        elif name == "RemoveHydrogens":
            stage = "Preparation"
            params = Chem.AddHsParameters()
            params.explicitOnly = False
            params.addCoords = False
            params.addResidueInfo = False
            params.skipQueries = False
            mol = Chem.AddHs(mol, params)
            stage = "Operation"
            params = Chem.RemoveHsParameters()
            for key in ("removeDegreeZero", "removeHigherDegrees", "removeOnlyHNeighbors",
                        "removeIsotopes", "removeAndTrackIsotopes", "removeDummyNeighbors",
                        "removeDefiningBondStereo", "removeWithQuery", "updateExplicitCount",
                        "removeHydrides", "removeNontetrahedralNeighbors"):
                setattr(params, key, False)
            for key in ("removeWithWedgedBond", "removeMapped", "removeInSGroups",
                        "showWarnings", "removeNonimplicit"):
                setattr(params, key, True)
            mol = Chem.RemoveHs(mol, params, sanitize=options["sanitize"])
        elif name == "Coordinates2dDefault":
            rdDepictor.SetPreferCoordGen(False)
            rdDepictor.Compute2DCoords(mol, canonOrient=False, clearConfs=True,
                coordMap={}, nFlipsPerSample=0, nSample=0, sampleSeed=0,
                permuteDeg4Nodes=False, forceRDKit=False, useRingTemplates=False)
            conf = mol.GetConformer()
            return {"Coordinates2d": {"topology": topology(mol), "xy_bits": [
                [struct.unpack(">Q", struct.pack(">d", coord))[0] for coord in
                 (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)]
                for i in range(mol.GetNumAtoms())]}}
        elif name != "SmilesRead":
            raise RuntimeError(f"unregistered molecular profile: {name}")
        return {"Topology": topology(mol)}
    except (ValueError, RuntimeError) as error:
        return {"Error": {"stage": stage, "detail": f"{type(error).__name__}: {error}"}}


results = []
for envelope in json.load(sys.stdin):
    if "Molecular" in envelope:
        results.append({"input": envelope, "output": {"Molecular": molecular(envelope["Molecular"])}})
        continue
    row = envelope["Fingerprint"]
    ctor = {"U32": DataStructs.UIntSparseIntVect,
            "U64": DataStructs.ULongSparseIntVect}[row["width"]]
    case = row["case"]
    def build(entries):
        value = ctor(case["length"])
        for key, count in entries:
            if count == 0:
                value[key] = 1
        value -= 1
        for key, count in entries:
            if count != 0:
                value[key] = count
        return value
    left, right = build(case["left"]), build(case["right"])
    before = (dict(left.GetNonzeroElements()), dict(right.GetNonzeroElements()))
    op = row["operation"]
    if op == "FuzzyAnd":
        result = left & right
    elif op == "FuzzyOr":
        result = left | right
    else:
        raise ValueError(op)
    assert before == (dict(left.GetNonzeroElements()), dict(right.GetNonzeroElements()))
    results.append({"input": envelope, "output": {"Fingerprint": {
        "length": result.GetLength(),
        "entries": sorted(result.GetNonzeroElements().items()),
    }}})
json.dump(results, sys.stdout, sort_keys=True)
