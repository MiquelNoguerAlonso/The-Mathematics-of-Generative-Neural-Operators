import GNO.FiniteCoupling
import Mathlib.Tactic

namespace GNO

/-- The Gram-matrix distortion for the two-coordinate nonorthogonal example. -/
theorem gram_distortion (c rho x y : ℝ) (hr : |c| ≤ rho) :
    (1-rho)*(x^2+y^2) ≤ x^2+y^2+2*c*x*y ∧
    x^2+y^2+2*c*x*y ≤ (1+rho)*(x^2+y^2) := by
  have hr0 : 0 ≤ rho := le_trans (abs_nonneg c) hr
  have hcross : |2*c*x*y| ≤ rho*(x^2+y^2) := by
    calc
      |2*c*x*y| = |c| *(2* |x| * |y|) := by simp only [abs_mul]; norm_num; ring
      _ ≤ rho*(2* |x| * |y|) := mul_le_mul_of_nonneg_right hr (by positivity)
      _ ≤ rho*(x^2+y^2) := by
        apply mul_le_mul_of_nonneg_left _ hr0
        nlinarith [sq_nonneg (|x|-|y|), sq_abs x, sq_abs y]
  rcases abs_le.mp hcross with ⟨hl, hu⟩
  constructor <;> nlinarith

/-- Transport the metric distortion to an actual finite coupling cost. -/
theorem finite_decoder_transfer {Ω : Type*} [Fintype Ω]
    (P : FiniteLaw Ω) (source target : Ω → ℝ) (a b : ℝ)
    (hlo : ∀ i, a*source i ≤ target i)
    (hhi : ∀ i, target i ≤ b*source i) :
    a*(∑ i, P.mass i*source i) ≤ ∑ i, P.mass i*target i ∧
    (∑ i, P.mass i*target i) ≤ b*(∑ i, P.mass i*source i) := by
  constructor
  · calc
      _ = ∑ i, P.mass i*(a*source i) := by simp only [Finset.mul_sum]; congr 1; funext i; ring
      _ ≤ _ := Finset.sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (hlo i) (P.nonneg i))
  · calc
      _ ≤ ∑ i, P.mass i*(b*source i) :=
        Finset.sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (hhi i) (P.nonneg i))
      _ = _ := by simp only [Finset.mul_sum]; congr 1; funext i; ring

/-- A finite law cannot approximate both point masses better than the
two-point lower bound on an unobserved feature cell of mass p. -/
theorem unseen_cell_two_point {Ω : Type*} [Fintype Ω]
    (P : FiniteLaw Ω) (z : Ω → ℝ) (p D : ℝ) (hp : 0 ≤ p) :
    p*D^2/4 ≤ max (p*∑ i, P.mass i*(z i)^2)
      (p*∑ i, P.mass i*(z i-D)^2) := by
  have hpoint : ∀ i, D^2/2 ≤ (z i)^2+(z i-D)^2 := by
    intro i; nlinarith [sq_nonneg (z i-D/2)]
  have hs := Finset.sum_le_sum (fun i (_ : i ∈ Finset.univ) =>
    mul_le_mul_of_nonneg_left (hpoint i) (P.nonneg i))
  simp only [mul_add, Finset.sum_add_distrib, ← Finset.sum_mul, P.total, one_mul] at hs
  have hm := mul_le_mul_of_nonneg_left hs hp
  have h1 := le_max_left (p*∑ i, P.mass i*(z i)^2) (p*∑ i, P.mass i*(z i-D)^2)
  have h2 := le_max_right (p*∑ i, P.mass i*(z i)^2) (p*∑ i, P.mass i*(z i-D)^2)
  linarith

end GNO
