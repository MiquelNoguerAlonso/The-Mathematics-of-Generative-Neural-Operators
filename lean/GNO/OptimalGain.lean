import GNO.Refinement
import Mathlib.Tactic

namespace GNO

/-- Largest eigenvalue of the real symmetric two by two matrix. -/
noncomputable def topEigenvalue (a b c : ℝ) : ℝ :=
  (a + c + Real.sqrt ((a-c)^2 + 4*b^2))/2

theorem topEigenvalue_diagonal (a b c : ℝ) :
    a ≤ topEigenvalue a b c ∧ c ≤ topEigenvalue a b c := by
  have hs := Real.sq_sqrt (show 0 ≤ (a-c)^2+4*b^2 by positivity)
  have hn := Real.sqrt_nonneg ((a-c)^2+4*b^2)
  have hdiff : |a-c| ≤ Real.sqrt ((a-c)^2+4*b^2) := by
    nlinarith [sq_abs (a-c), abs_nonneg (a-c), sq_nonneg b]
  unfold topEigenvalue
  rcases abs_le.mp hdiff with ⟨hlo, hhi⟩
  constructor <;> linarith

theorem topEigenvalue_characteristic (a b c : ℝ) :
    (topEigenvalue a b c-a)*(topEigenvalue a b c-c) = b^2 := by
  have hs := Real.sq_sqrt (show 0 ≤ (a-c)^2+4*b^2 by positivity)
  unfold topEigenvalue
  nlinarith

/-- The optimized quadratic bound, with no sign restriction on the vector. -/
theorem topEigenvalue_quadratic (a b c x y : ℝ) :
    a*x^2+2*b*x*y+c*y^2 ≤ topEigenvalue a b c*(x^2+y^2) := by
  let g := topEigenvalue a b c
  have hd := topEigenvalue_diagonal a b c
  have hc := topEigenvalue_characteristic a b c
  change a ≤ g ∧ c ≤ g at hd
  change (g-a)*(g-c)=b^2 at hc
  change a*x^2+2*b*x*y+c*y^2 ≤ g*(x^2+y^2)
  by_cases heq : g = a
  · have hb : b = 0 := by rw [heq] at hc; nlinarith [sq_nonneg b]
    rw [hb]
    nlinarith [mul_nonneg (sub_nonneg.mpr hd.2) (sq_nonneg y)]
  · have hp : 0 < g-a := by rcases hd with ⟨ha, _⟩; exact sub_pos.mpr (lt_of_le_of_ne ha (Ne.symm heq))
    have hid : (g-a)*(g*(x^2+y^2)-(a*x^2+2*b*x*y+c*y^2)) =
        ((g-a)*x-b*y)^2 := by
      nlinarith [congrArg (fun z : ℝ => z*y^2) hc]
    have hm : 0 ≤ (g-a)*(g*(x^2+y^2)-(a*x^2+2*b*x*y+c*y^2)) := by
      rw [hid]; exact sq_nonneg _
    exact sub_nonneg.mp ((mul_nonneg_iff_of_pos_left hp).mp hm)

/-- A nonzero vector attains the bound, including diagonal degeneracies. -/
theorem topEigenvalue_attained (a b c : ℝ) :
    ∃ x y : ℝ, 0 < x^2+y^2 ∧
      a*x^2+2*b*x*y+c*y^2 = topEigenvalue a b c*(x^2+y^2) := by
  let g := topEigenvalue a b c
  have hd := topEigenvalue_diagonal a b c
  have hc := topEigenvalue_characteristic a b c
  change a ≤ g ∧ c ≤ g at hd
  change (g-a)*(g-c)=b^2 at hc
  change ∃ x y : ℝ, 0 < x^2+y^2 ∧ a*x^2+2*b*x*y+c*y^2=g*(x^2+y^2)
  by_cases heq : g = a
  · exact ⟨1, 0, by norm_num, by simp [heq]⟩
  · refine ⟨b, g-a, ?_, ?_⟩
    · have hs := sq_pos_of_ne_zero (sub_ne_zero.mpr heq)
      nlinarith [sq_nonneg b]
    · nlinarith [congrArg (fun z : ℝ => z*(g-a)) hc]

/-- Exact matrix bound used at one stage of the aggregate-energy recursion. -/
theorem optimal_energy_step (g L x e : ℝ) (hg : 0 ≤ g) :
    g*x^2+(e+L*Real.sqrt g*x)^2 ≤
      topEigenvalue ((1+L^2)*g) (L*Real.sqrt g) 1 * (x^2+e^2) := by
  have hs := Real.sq_sqrt hg
  have hb := topEigenvalue_quadratic ((1+L^2)*g) (L*Real.sqrt g) 1 x e
  nlinarith [congrArg (fun z : ℝ => z*(L*x)^2) hs]

end GNO
