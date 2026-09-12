# The Mathematics of Generative Neural Operators

**Manuscript:** [Read the complete PDF](Generative_Neural_Operators.pdf)

**Author:** Miquel Noguer i Alonso

**DOI:** [10.5281/zenodo.22714854](https://doi.org/10.5281/zenodo.22714854)

The manuscript develops resolution-autonomy obstructions, a Gaussian transport-cost premium, conditional refinement generators, a sharp Wasserstein error recurrence, an infinite-dimensional existence certificate, and observation and physical-time limits.

## Build the paper

Main file: `Generative_Neural_Operators.tex`.

The source uses `natbib` with author-year citations and `\bibliographystyle{plainnat}`. All seven figures are PNG files. The numerical tables and figures are included, so compiling does not require running Python.

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error Generative_Neural_Operators.tex
```

Alternatively run `pdflatex`, `bibtex`, then `pdflatex` twice. For Overleaf, upload the complete source ZIP and select `Generative_Neural_Operators.tex` as the main document and pdfLaTeX as the compiler.

## Reproduce the numerical illustrations

```bash
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiments/reproduce.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiments/optimal_gain.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiments/learn_operator.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiments/observed_audit.py
```

The first script regenerates four figures, two LaTeX tables, and `results/metrics.json`, which records its environment and numerical results. The synthetic random-field experiment uses seed 20260911, 60,000 independent realizations, and a 256-mode finite reference. Gaussian matrix ODEs use DOP853 with rtol 1e-11 and atol 1e-13.

The first script contains the original explicit-kernel experiments. The optimal-gain script computes the exact energy amplification and reconstructs attaining error vectors. The learned-operator script trains three preselected 37-parameter shared spectral networks and a context ablation on noisy field coefficients, then performs an independent conditional audit. No real-data benchmark superiority is claimed. Nonlinear coupling RMS values are distinguished from optimal Wasserstein distances and from deterministic theorem bounds. All examples use synthetic data and require no credentials or external datasets.

## Files

- `Generative_Neural_Operators.tex`: main LaTeX document.
- `sections/`: complete manuscript text and proofs.
- `references.bib`: verified bibliography using plainnat.
- `figures/`: seven publication figures in PNG format.
- `tables/`: generated LaTeX tables.
- `experiments/reproduce.py`: self-contained numerical reproduction and formula checks.
- `experiments/optimal_gain.py`: exact energy gains and attaining errors.
- `experiments/learn_operator.py`: training, ablation, and independent validation.
- `models/`: seven learned parameter files in portable JSON.
- `results/metrics.json`: full results, sampling errors, and package versions.
- `results/learned_operator.json` and `results/optimal_gain.json`: new experimental records.
- `results/revision_validation.json`: executed reproduction and proof-audit summary.
- `Generative_Neural_Operators.pdf`: compiled manuscript.

The results concern the stated mathematical model classes and population errors. The Gaussian audit uses a known conditional mean. The bounded-detail audit uses observed fields and independent aggregation, under specified support, sufficient-context, regularity, and tail assumptions. Optimization guarantees and performance on market or PDE datasets remain outside its scope.

## Lean verification

The lean directory is an independently buildable Lean 4.19.0 project pinned to mathlib v4.19.0. Its 35 theorem declarations pass compilation and an axiom audit. The project proves the finite conditional regression identity, real Hilbert orthogonality identities, the finite refinement comparison and exponential bound, properties of the stated Gaussian midpoint formula, a metric two-point bound, and a fiber-factorization theorem.

~~~bash
cd lean
lake exe cache get
python3 verify.py
~~~

The coverage map, build log, axiom log, exact dependency manifest, and machine-readable source-hash report are included. Appendix A explains the precise formalization boundary. In particular, the probability-kernel coupling construction and the infinite-dimensional existence results are written proofs, not fully mechanized results. There are no admitted proofs or custom mathematical axioms in the supplied Lean code.

## Learned model and validation protocol

The networks use 8,192 independent training fields through 32 modes, 1,800 Adam steps, and projected global sensitivity constraints. An independent audit uses 65,536 fields through 256 modes. The simultaneous 99% statement covers all three full-context models and 255 detail bands; the separate context ablation is not part of that joint confidence statement. Training seeds and hyperparameters are fixed before the audit.

CPU float64 and deterministic single-threaded PyTorch are used. All weights are saved as portable JSON in models/, without executable pickles. Results are in results/learned_operator.json and results/optimal_gain.json. Coupling costs, plug-in recurrence estimates, confidence bounds, and unresolved tails are distinct quantities. The known synthetic conditional mean is used only by the audit, not as a training label.

The new Lean module derives a finite weighted joint-coupling step from its actual probability law. The isolated-seed theorem checks the finite product underlying the sharp robustness threshold. The exact two-dimensional spectral upper bound, characteristic identity, attainment, and energy specialization are checked in Lean. Optimal measurable coupling, hierarchy-wide attainment, concentration, and infinite limits remain written proofs.

## Observed-sample audit and observation geometry

`experiments/observed_audit.py` trains three 24-parameter models on 16,384 noisy fields through eight modes. Independent calibration and aggregation splits audit 31 later modes using three predeclared budgets. The joint 99% event covers every candidate, cutoff, and budget. The audit does not call the target conditional mean. It uses exact one-dimensional empirical-to-uniform transport, a bounded-support DKW certificate, and Hoeffding aggregation. A separate truth diagnostic evaluates a specified coupling.

`results/observed_audit.json` records both the generic prefix recurrence and the tighter dependency-sensitive certificate, all confidence inputs and diagnostic values. Three additional model JSON files, one table, and two figures are generated. Riesz-coordinate and bi-Lipschitz decoder results identify how the physical metric changes. Lean additionally checks triangular dependency propagation, two-coordinate Gram bounds, finite coupling metric transfer, and a finite unseen-cell lower bound. Statistical coverage and floating-point arithmetic are distinct claims.

## Standing academic LaTeX format

Use the author's established academic format for this and future papers: 11pt article; 1.08-inch margins; 1.08 line spacing; 1.25em paragraph indentation and no paragraph skip; small captions with bold labels; concealed link styling; natbib/plainnat; the Miquel Noguer i Alonso / Artificial Intelligence Finance Institute (AIFI) author block; and an actual `\today` date. The front matter includes a clickable table of contents covering sections and subsections, the appendix, and references. The contents and main text begin on fresh pages. Running headers, footers, and printed page numbers remain omitted. DOI links are separate title-page identifiers. `Academic_LaTeX_Format.md` records the reusable specification.
