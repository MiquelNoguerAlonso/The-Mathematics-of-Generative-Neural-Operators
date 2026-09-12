import Mathlib.Tactic

namespace GNO

section FiniteRegression
variable {I : Type*} [Fintype I]

/-- The weighted conditional mean on one observation fiber. -/
noncomputable def fiberMean (w v : I → ℝ) : ℝ :=
  (∑ i, w i * v i) / (∑ i, w i)

theorem fiber_centered (w v : I → ℝ) (hm : (∑ i, w i) ≠ 0) :
    ∑ i, w i * (v i - fiberMean w v) = 0 := by
  simp_rw [mul_sub]
  rw [Finset.sum_sub_distrib, ← Finset.sum_mul]
  unfold fiberMean
  field_simp

/-- Exact loss decomposition for every finite fiber, not a sampled numerical
check. Nonnegative weights are only needed for the subsequent inequality. -/
theorem fiber_regression_identity (w v : I → ℝ) (b : ℝ)
    (hm : (∑ i, w i) ≠ 0) :
    (∑ i, w i * (v i-b)^2) =
      (∑ i, w i * (v i-fiberMean w v)^2) +
      (∑ i, w i) * (fiberMean w v-b)^2 := by
  have hpoint : ∀ i, w i * (v i-b)^2 =
      w i * (v i-fiberMean w v)^2 +
      w i * (fiberMean w v-b)^2 +
      (2*(fiberMean w v-b)) * (w i*(v i-fiberMean w v)) := by
    intro i
    ring
  simp_rw [hpoint]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib,
      ← Finset.sum_mul, ← Finset.mul_sum, fiber_centered w v hm]
  ring

theorem fiber_mean_minimizes (w v : I → ℝ) (b : ℝ)
    (hw : ∀ i, 0 ≤ w i) (hm : (∑ i, w i) ≠ 0) :
    (∑ i, w i * (v i-fiberMean w v)^2) ≤
      (∑ i, w i * (v i-b)^2) := by
  rw [fiber_regression_identity w v b hm]
  have hsum : 0 ≤ ∑ i, w i := Finset.sum_nonneg (fun i _ => hw i)
  have hp : 0 ≤ (∑ i, w i) * (fiberMean w v-b)^2 := by positivity
  linarith

/-- Conditional expectation loss decomposition on a finite coarse/fine
product space. Fibers must have positive mass; null fibers may be omitted. -/
theorem finite_conditional_regression {C : Type*} [Fintype C]
    (w v : C → I → ℝ) (b : C → ℝ)
    (hm : ∀ c, (∑ i, w c i) ≠ 0) :
    (∑ c, ∑ i, w c i * (v c i-b c)^2) =
      (∑ c, ∑ i, w c i * (v c i-fiberMean (w c) (v c))^2) +
      ∑ c, (∑ i, w c i) * (fiberMean (w c) (v c)-b c)^2 := by
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro c _
  exact fiber_regression_identity (w c) (v c) (b c) (hm c)

end FiniteRegression
end GNO
