# Generative Neural Operators

**Miquel Noguer i Alonso**  
Artificial Intelligence Finance Institute (AIFI)

Two mathematical papers on generative neural operators, with complete LaTeX sources, compiled PDFs, Lean 4 projects, reproducible synthetic experiments, trained parameters, and validation reports.

| Paper | Manuscript | DOI | Sources and reproduction |
|---|---|---|---|
| I. The Mathematics of Generative Neural Operators | [PDF, 49 pages](foundations/Generative_Neural_Operators.pdf) | [10.5281/zenodo.22714854](https://doi.org/10.5281/zenodo.22714854) | [Foundations README](foundations/README.md) |
| II. Generative Neural Operators in Finance | [PDF, 42 pages](finance/Generative_Neural_Operators_Finance.pdf) | [10.5281/zenodo.22714857](https://doi.org/10.5281/zenodo.22714857) | [Finance README](finance/README.md) |

Both papers include linked tables of contents, author-year citations with `plainnat`, and the author's academic LaTeX format. The source uses `\today`; subsequent builds display their compilation date.

## Paper I: mathematical foundations

The foundations paper studies when generative laws remain consistent across observation resolutions and why a projected state may fail to admit autonomous dynamics. Conditional refinement preserves generated coarse coordinates. A sharp orthogonal error recurrence separates conditional fitting error from unresolved energy, while an exact energy gain characterizes class-uniform robustness.

The paper also treats nonorthogonal observations, invertible nonlinear representations, architecture-dependent error propagation, and an independent observed-sample audit under explicit support, regularity, sufficient-context, and tail assumptions.

## Paper II: finance

The finance paper connects latent field approximation, posterior inference, and event-intensity errors to observable market paths and execution decisions. It studies hidden liquidity, reserve conservation, filtering, policy comparison, and martingale pricing.

For a finite observable control model, recorded counts and exposures produce simultaneous rate intervals. Robust dynamic programming turns those intervals into guarantees for a policy selected from the same data. Continuation values quantify economically relevant rate errors, numerical residuals account for time discretization, and an early-exercise example isolates the importance of information timing.

## Build the PDFs

Each paper is self-contained. Its figures, tables, and bibliography are already included; building the PDF does not require running Python or Lean.

From the repository root:

```bash
cd foundations
latexmk -pdf -interaction=nonstopmode -halt-on-error Generative_Neural_Operators.tex
cd ../finance
latexmk -pdf -interaction=nonstopmode -halt-on-error Generative_Neural_Operators_Finance.tex
```

Use pdfLaTeX and BibTeX. For Overleaf, upload the contents of the relevant paper folder and select its main `.tex` file.

## Run the Lean projects

Both projects pin Lean 4.19.0 and their mathlib dependency. With the Lean toolchain manager installed, run:

```bash
cd foundations/lean
lake exe cache get
python3 verify.py
cd ../../finance/lean
lake exe cache get
python3 verify.py
```

The supplied audits report **35 declarations for foundations and 40 for finance**, with no admitted proofs or custom mathematical axioms. Ten refinement declarations are shared between the standalone projects; the total is not a count of distinct new research results.

The projects mechanize selected finite probabilistic, algebraic, and geometric results. General disintegration, infinite-dimensional limits, concentration results, continuous-time stochastic control, and bicausal stopping remain written proofs. Read the declaration-level [foundations coverage map](foundations/lean/COVERAGE.md) and [finance coverage map](finance/lean/COVERAGE.md) for the exact boundaries.

## Reproduce the experiments

Install the dependencies from the relevant paper's `requirements.txt`, then follow its README. The scripts and saved parameters use synthetic data and require no credentials or external market datasets.

| Folder | Programs |
|---|---|
| [Foundations experiments](foundations/experiments/) | Explicit-kernel illustrations, exact energy gains, shared spectral network training, and observed-sample conditional audits |
| [Finance experiments](finance/experiments/) | Filtering and accounting examples, a learned reserve-rate model, and event-count-based inference with robust execution |

Each paper includes portable model parameters in `models/`, numerical reports in `results/`, PNG figures, and generated LaTeX tables. The supplied validation records document exact reruns in the recorded environment. Floating-point calculations and empirical diagnostics are separate from exact Lean proofs and statistical coverage statements. The experiments do not establish exchange-data calibration or trading profitability.

## Review and verification records

- [Technical assessment of both papers](review/Technical_Assessment.md): an internal technical review, with the reviewer's role and remaining scientific limitations disclosed.
- [Foundations document audit](foundations/results/document_audit.json) and [finance document audit](finance/results/document_audit.json): PDF, citation, contents-link, and layout checks.
- [Foundations reproduction record](foundations/results/revision_validation.json) and [finance reproduction record](finance/results/revision_validation.json): executed numerical checks and source hashes.
- [Foundations Lean audit](foundations/lean/verification.json) and [finance Lean audit](finance/lean/verification.json): checked declarations, axiom dependencies, and source hashes.

`MANIFEST_SHA256.txt` inventories the repository snapshot. Each paper also retains its own package inventory and complete build instructions.
