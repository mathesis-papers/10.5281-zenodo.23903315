# Shannon branch-and-bound counterexample

This note refutes the Shannon conjecture following Eq. (69), Section V.C,
of Tarabunga et al., arXiv:2411.05766v1 (Eq. (69) on p. 9, the explicit
Shannon conjecture on p. 10). The signed witness is encoded
in the linked Trureturing declaration. The note also supplies written proofs
for nonnegative amplitudes, strictly positive existence, every dimension
at least three, a Rényi-index interval, and a collision-entropy replacement
bound for arbitrary affine supports.
It also gives the explicit phase conversion for the real-amplitude reduction
in Appendix C.2: integer-sum linear phases require a corresponding quadratic
adjustment. The global reduction to real bases remains valid.
The note determines the all-index minimum for the witness at fixed `R=I3`,
proves extensive growth of the Shannon defect, strengthens the collision
bound by conditioning on affine cosets, and proves validity of the original
bound on coset-flat amplitudes. In particular, that bound is valid for the
posing source's Figure 1 product-state families; this is not a statement
about the numerical implementation.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
python3 verify.py --output verification.json
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is contained
in `note.tex`; no BibTeX pass or separate bibliography file is needed.
Python 3.10 or later suffices for the standard-library finite checks.
No other local source files are required to build the PDF.

## Mathematical scope

- The pinned Lean declaration `D5.S3.Quantum.Magic.StabilizerAsymmetryShannonBoundRefutation.result`
  proves the negation of the full-support Shannon inequality, witnessed by
  the signed three-qubit state. Its source is linked in the PDF and metadata.
- The additional corollaries and arbitrary-support collision bound have
  complete written proofs in the note.
- `verify.py` computes the transforms exactly, checks all 64 full-support
  quadratic/linear phase choices for the nonnegative witness, and certifies
  the strictly positive example by rational logarithm intervals.
- The fixed-support minimum is `ln(4)`. The calculation does not minimize
  over other support matrices and does not determine the full BMSA.
- The all-dimensional extension is a written tensor-product proof; the
  script also checks padding in dimensions 3 through 8.

## Contents and license

`note.tex`, `note.pdf`, `metadata.json`, `verify.py` and `verification.json`
form the scientific package. `source-v1.zip` contains these text sources
and this README, with a SHA-256 member manifest and without the PDF.

Creators: Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).

License: Creative Commons Attribution 4.0 International (CC BY 4.0),
identifier `cc-by-4.0`.
