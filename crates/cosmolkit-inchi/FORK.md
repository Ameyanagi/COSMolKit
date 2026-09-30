# ReShiki compatibility branch

This branch starts from upstream COSMolKit 0.3.0, commit
`d892ec3507c5b568c5ed5d86ae44e466f7d03855`, and retains the upstream MIT license.
It carries three corrections used by ReShiki's Rust InChI helper:

- Restore implicit hydrogen addition for Tl through Ra from the official
  InChI 1.07.5 `ElData[].bSkipAddingH` table. The adjacent Hg and Ac ranges
  still suppress automatic addition.
- Preserve bond insertion order when enumerating perchlorate cleanup
  matches, as RDKit's substructure matcher does. Sorting neighbors by atom
  index can change which equivalent substructure is rewritten.

- Emit explicit unknown-parity stereo records for `BondStereo::Any`, including
  when stereo atom annotations are absent, matching RDKit 2026.03.6 commit
  `0e0d85f4ca34aeae15dfc0f7cf5503bdb0a8e985`. Its BSD notice is retained in
  `LICENSE-RDKit`.

Focused regressions live beside each correction. ReShiki also compares
independently captured molecular inputs with the official InChI 1.07.5
kernel, commit `11a87982bb518f57ac013f0b258c283655e1ea1d`.

The native reference submodules are configured with `update = none` so Cargo
Git dependency resolution does not fetch unrelated C/C++ source trees. They
are unnecessary to build this Rust crate. Developers running native oracle
tools can fetch them explicitly:

```sh
git submodule update --init --checkout third_party/rdkit third_party/InChI
```

Run the focused tests with:

```sh
cargo test -p cosmolkit-inchi --lib inchi_compatibility
cargo test -p cosmolkit-inchi --lib source_port__util__
cargo test -p cosmolkit-inchi --lib source_port__inchi__rcleanup__line_1694
```
