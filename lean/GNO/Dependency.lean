import Mathlib.Tactic

namespace GNO

noncomputable def dependencyRadius (eps : ℕ → ℝ) (B : ℕ → ℕ → ℝ)
    (n : ℕ) : ℝ := eps n + ∑ i : Fin n, B n i * dependencyRadius eps B i
termination_by n

/-- Triangular propagation preserves a coordinatewise conditional-error bound.
Its probabilistic coupling hypotheses are supplied by the manuscript. -/
theorem dependency_comparison (eps D : ℕ → ℝ) (B : ℕ → ℕ → ℝ)
    (hB : ∀ n i, 0 ≤ B n i)
    (hstep : ∀ n, D n ≤ eps n+∑ i : Fin n, B n i*D i) :
    ∀ n, D n ≤ dependencyRadius eps B n := by
  intro n
  induction n using Nat.strong_induction_on with
  | h n ih =>
      rw [dependencyRadius]
      apply le_trans (hstep n)
      apply add_le_add_left
      exact Finset.sum_le_sum (fun i _ =>
        mul_le_mul_of_nonneg_left (ih i i.isLt) (hB n i))

theorem dependency_radius_nonneg (eps : ℕ → ℝ) (B : ℕ → ℕ → ℝ)
    (he : ∀ n, 0 ≤ eps n) (hB : ∀ n i, 0 ≤ B n i) :
    ∀ n, 0 ≤ dependencyRadius eps B n := by
  apply dependency_comparison eps (fun _ => 0) B hB
  intro n
  simpa using he n

end GNO
