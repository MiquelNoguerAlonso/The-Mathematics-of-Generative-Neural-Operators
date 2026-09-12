# Formalization coverage

**Verified theorem declarations: 35.** This includes supporting lemmas. The refinement module appears in both standalone projects.

The Lean build and per-theorem axiom audit pass. There are no admitted proofs or custom mathematical axioms. This is a partial formalization of the manuscript, with the exact scope below.

## Manuscript correspondence

| Manuscript component | Formal coverage and boundary |
| --- | --- |
| Exact autonomy defect and hierarchy | Finite scalar conditional-regression analogue and real Hilbert orthogonality; general conditional expectation and time integrals are written proofs. |
| Gaussian covariance-commutation theorem | Only the stated two-dimensional midpoint formula properties are checked; covariance regression and matrix commutation are written proofs. |
| Coherent conditional representation | Written proof only: disintegration, kernel randomization, projective sampling, and Hilbert limits are not mechanized. |
| Gaussian triangular transport premium | Written proof only; the Gaussian optimal transport and kinetic-action arguments are not mechanized. |
| Orthogonal refinement certificate | The finite real-sequence comparison and exponential bound are checked. Wasserstein couplings, exact tail transport identity, and infinite limit are written proofs. |
| Sharpness and convergence corollaries | The recurrence itself is checked; the Dirac-kernel construction and infinite law convergence are written proofs. |
| Flow-matching loss-to-law result | Finite regression analogue only; ODE stability, continuous-time integration, and numerical-solver bounds are written proofs. |
| Minimax unresolved completion | The abstract two-point metric step is checked. The indistinguishable probability-law construction is a written proof. |
| Autonomous physical-time kernel criterion | Pure fiber-factorization equivalence checked; measurable quotient and randomized kernel constructions are written proofs. |
| Convex constraints and projected laws | Written proof only; metric projection and Wasserstein pushforward statements are not mechanized. |

## Every checked declaration

| Declaration | Source | Exact checked claim |
| --- | --- | --- |
| GNO.orthogonal_add_square | GNO/Geometry.lean | Pythagorean identity in an arbitrary real inner-product space. |
| GNO.orthogonal_three_square | GNO/Geometry.lean | Three mutually orthogonal vectors have additive squared norms. |
| GNO.midpointDefect_nonneg | GNO/Geometry.lean | Nonnegativity of the defined two-dimensional Gaussian midpoint formula for rho squared below one. |
| GNO.midpointDefect_pos_iff | GNO/Geometry.lean | The defined midpoint formula is positive exactly for nonzero correlation, under rho squared below one. |
| GNO.midpointDefect_four_fifths | GNO/Geometry.lean | Exact value 8/21 at correlation 4/5. |
| GNO.two_point_lower_bound | GNO/Information.lean | Universal metric two-target lower bound against every common estimate. |
| GNO.factors_iff_fiber_constant | GNO/Information.lean | Factorization through a surjective observation map exactly when the map is constant on fibers; no measurability assertion. |
| GNO.radius_nonneg | GNO/Refinement.lean | The recursively defined real radius is nonnegative at every depth. |
| GNO.radius_sq_succ | GNO/Refinement.lean | The square-root definition satisfies the exact squared recurrence. |
| GNO.radius_monotone | GNO/Refinement.lean | The radius is nondecreasing with the number of retained bands. |
| GNO.square_step_monotone | GNO/Refinement.lean | The squared update is monotone on nonnegative errors for nonnegative fitting error and sensitivity. |
| GNO.error_le_radius | GNO/Refinement.lean | Every nonnegative discrepancy sequence satisfying the one-step squared bound lies below the radius. |
| GNO.radius_step_bound | GNO/Refinement.lean | The nonlinear squared step is bounded by its linear majorant step. |
| GNO.radius_sq_le_linearBound | GNO/Refinement.lean | The radius squared lies below the recursively defined linear majorant. |
| GNO.context_independent_exact | GNO/Refinement.lean | Zero context sensitivity gives the exact sum of squared fitting errors. |
| GNO.exponential_certificate | GNO/Refinement.lean | Universal finite exponential bound on the radius squared. |
| GNO.add_tail_certificate | GNO/Refinement.lean | Algebraic consequence of a supplied orthogonal square decomposition and projected square bound. |
| GNO.fiber_centered | GNO/Regression.lean | Weighted deviations from the finite fiber mean sum to zero. |
| GNO.fiber_regression_identity | GNO/Regression.lean | Exact weighted squared-loss decomposition on one nonzero-mass finite fiber. |
| GNO.fiber_mean_minimizes | GNO/Regression.lean | The weighted fiber mean minimizes squared loss for nonnegative weights. |
| GNO.finite_conditional_regression | GNO/Regression.lean | Exact scalar conditional-expectation regression identity on arbitrary finite coarse/fine product spaces. |

## Spectral optimization, dependency, and observation extensions

| Declaration | Checked statement | Boundary |
|---|---|---|
| GNO.topEigenvalue_diagonal | Largest symmetric 2x2 eigenvalue dominates both diagonal entries | Exact real formula |
| GNO.topEigenvalue_characteristic | Characteristic polynomial identity | Exact real formula |
| GNO.topEigenvalue_quadratic | Universal quadratic-form upper bound | Arbitrary real vector |
| GNO.topEigenvalue_attained | A nonzero vector attains the bound, including degeneracy | Hierarchy-wide kernel attainment remains analytic |
| GNO.optimal_energy_step | Spectral bound specialized to one aggregate-energy update | Coupling and complete optimal-gain induction remain analytic |
| GNO.dependency_comparison | Well-founded triangular propagation dominates any sequence satisfying the local inequalities | Probabilistic hypotheses supplied separately |
| GNO.dependency_radius_nonneg | Nonnegative input and dependency coefficients give nonnegative radii | Arbitrary finite cutoff |
| GNO.gram_distortion | Both exact two-coordinate Gram inequalities | General Riesz operators remain analytic |
| GNO.finite_decoder_transfer | Pointwise metric bounds transfer to a finite normalized coupling cost | Infima over general couplings remain analytic |
| GNO.unseen_cell_two_point | Finite output law has the stated two-point squared-loss lower bound | Sampling probability and confidence consequence remain analytic |

The observed-sample confidence theorem, DKW--Massart inequality, conditional cell laws, independent aggregation, and bicausal/measure-theoretic constructions are not mechanized by these declarations. The Python audit does not constitute a formal proof of its floating-point output.

## Trust and reproduction

The statements use exact real numbers and finite types, not floating-point experiments. Standard library foundations are reported for every theorem in axioms.log. verification.json records their axiom sets and the SHA-256 of each checked source file and project configuration.

Hypotheses appear in the Lean declarations as hypotheses; they are not disguised as proofs. Probabilistic and analytic constructions excluded above remain explicitly written proofs in the manuscript. A theorem declaration in this project is never an assertion that every result in its associated paper section has been mechanized.

## Added finite probability and sharpness interfaces

| Declaration | Checked statement | Analytic boundary |
|---|---|---|
| GNO.weighted_cross_bound | Weighted finite Cauchy--Schwarz with second-moment bounds | No general measure or optimal transport construction |
| GNO.weighted_minkowski_square | Squared finite weighted Minkowski bound | No disintegration |
| GNO.finite_refinement_step | Actual finite joint coupling cost under fitting and pointwise context bounds | Optimal couplings and general kernel measurability remain written proofs |
| GNO.isolated_seed_product | Exact finite product when fitting error occurs only in the first band | Infinite product divergence and the optimal eigenvalue gain remain written proofs |

The shared spectral architecture and empirical Bernstein validation corollary are analytic results in the manuscript; they are not mechanized here.
