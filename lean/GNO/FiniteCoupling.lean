import GNO.Geometry
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Tactic

namespace GNO

/-- An actual finite probability law, not an assumed error sequence. -/
structure FiniteLaw (Ω : Type*) [Fintype Ω] where
  mass : Ω → ℝ
  nonneg : ∀ i, 0 ≤ mass i
  total : ∑ i, mass i = 1

variable {Ω : Type*} [Fintype Ω]

theorem weighted_cross_bound (w a b : Ω → ℝ) (hw : ∀ i, 0 ≤ w i)
    {E F : ℝ} (hE : 0 ≤ E) (hF : 0 ≤ F)
    (ha : ∑ i, w i * (a i)^2 ≤ E^2)
    (hb : ∑ i, w i * (b i)^2 ≤ F^2) :
    (∑ i, w i * a i * b i) ≤ E*F := by
  have hca := Finset.sum_sq_le_sum_mul_sum_of_sq_eq_mul
    (s := Finset.univ)
    (r := fun i => w i * a i * b i)
    (f := fun i => w i * (a i)^2)
    (g := fun i => w i * (b i)^2)
    (by intro i _; exact mul_nonneg (hw i) (sq_nonneg _))
    (by intro i _; exact mul_nonneg (hw i) (sq_nonneg _))
    (by intro i _; ring)
  have hag : 0 ≤ ∑ i, w i * (a i)^2 :=
    Finset.sum_nonneg (by intro i _; exact mul_nonneg (hw i) (sq_nonneg _))
  have hbg : 0 ≤ ∑ i, w i * (b i)^2 :=
    Finset.sum_nonneg (by intro i _; exact mul_nonneg (hw i) (sq_nonneg _))
  have hp := mul_le_mul ha hb hbg (sq_nonneg E)
  have hef : 0 ≤ E*F := mul_nonneg hE hF
  nlinarith [sq_nonneg ((∑ i, w i * a i * b i) - E*F)]

theorem weighted_minkowski_square (w a b : Ω → ℝ) (hw : ∀ i, 0 ≤ w i)
    {E F : ℝ} (hE : 0 ≤ E) (hF : 0 ≤ F)
    (ha : ∑ i, w i * (a i)^2 ≤ E^2)
    (hb : ∑ i, w i * (b i)^2 ≤ F^2) :
    ∑ i, w i * (a i + b i)^2 ≤ (E+F)^2 := by
  have hc := weighted_cross_bound w a b hw hE hF ha hb
  have heq : (∑ i, w i * (a i+b i)^2) =
      (∑ i, w i*(a i)^2) + 2*(∑ i, w i*a i*b i) +
      (∑ i, w i*(b i)^2) := by
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [heq]
  nlinarith

/-- The finite probabilistic one-step interface for the orthogonal certificate.
X,Y are coupled old contexts and A,B,C are the true, intermediate, and learned
scalar details. The conclusion bounds the actual joint squared coupling cost.
Optimality and measurable construction of such couplings remain analytic. -/
theorem finite_refinement_step {H : Type*} [NormedAddCommGroup H]
    (P : FiniteLaw Ω) (X Y : Ω → H) (A B C : Ω → ℝ)
    {eps L D : ℝ} (he : 0 ≤ eps) (hL : 0 ≤ L) (hD : 0 ≤ D)
    (hold : ∑ i, P.mass i * ‖X i-Y i‖^2 ≤ D^2)
    (hfit : ∑ i, P.mass i * (A i-B i)^2 ≤ eps^2)
    (hcontext : ∀ i, |B i-C i| ≤ L*‖X i-Y i‖) :
    ∑ i, P.mass i * (‖X i-Y i‖^2+(A i-C i)^2) ≤
      D^2+(eps+L*D)^2 := by
  have hbpoint : ∀ i, (B i-C i)^2 ≤ L^2*‖X i-Y i‖^2 := by
    intro i
    have hn := norm_nonneg (X i-Y i)
    have ha := abs_nonneg (B i-C i)
    have hs := sq_abs (B i-C i)
    nlinarith [hcontext i, sq_nonneg (L*‖X i-Y i‖-|B i-C i|)]
  have hb : (∑ i, P.mass i*(B i-C i)^2) ≤ (L*D)^2 := by
    calc
      _ ≤ ∑ i, P.mass i*(L^2*‖X i-Y i‖^2) := by
        apply Finset.sum_le_sum
        intro i _
        exact mul_le_mul_of_nonneg_left (hbpoint i) (P.nonneg i)
      _ = L^2*(∑ i, P.mass i*‖X i-Y i‖^2) := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro i _
        ring
      _ ≤ L^2*D^2 := mul_le_mul_of_nonneg_left hold (sq_nonneg L)
      _ = (L*D)^2 := by ring
  have ht := weighted_minkowski_square P.mass
    (fun i => A i-B i) (fun i => B i-C i) P.nonneg
    he (mul_nonneg hL hD) hfit hb
  have hsimp : ∀ i, A i-B i+(B i-C i) = A i-C i := by
    intro i
    ring
  simp only at ht
  simp_rw [hsimp] at ht
  simp_rw [mul_add, Finset.sum_add_distrib]
  exact add_le_add hold ht

end GNO
