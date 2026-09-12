import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

namespace GNO

section Hilbert
variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

theorem orthogonal_add_square (x y : H)
    (h : inner (𝕜 := ℝ) x y = 0) :
    ‖x+y‖^2 = ‖x‖^2 + ‖y‖^2 := by
  rw [norm_add_sq_real, h]
  ring

theorem orthogonal_three_square (a b c : H)
    (hab : inner (𝕜 := ℝ) a b = 0)
    (hac : inner (𝕜 := ℝ) a c = 0)
    (hbc : inner (𝕜 := ℝ) b c = 0) :
    ‖a+b+c‖^2 = ‖a‖^2 + ‖b‖^2 + ‖c‖^2 := by
  have hsum : inner (𝕜 := ℝ) (a+b) c = 0 := by
    rw [inner_add_left, hac, hbc]
    ring
  rw [orthogonal_add_square (a+b) c hsum,
      orthogonal_add_square a b hab]

end Hilbert

noncomputable def midpointDefect (rho : ℝ) : ℝ :=
  2*rho^2/(4-rho^2)

theorem midpointDefect_nonneg {rho : ℝ} (h : rho^2 < 1) :
    0 ≤ midpointDefect rho := by
  unfold midpointDefect
  exact div_nonneg (by positivity) (by linarith)

theorem midpointDefect_pos_iff {rho : ℝ} (h : rho^2 < 1) :
    0 < midpointDefect rho ↔ rho ≠ 0 := by
  have hd : 0 < 4-rho^2 := by linarith
  unfold midpointDefect
  rw [div_pos_iff_of_pos_right hd]
  constructor
  · intro hp hz
    simp [hz] at hp
  · intro hn
    have hs : 0 < rho^2 := sq_pos_of_ne_zero hn
    positivity

theorem midpointDefect_four_fifths :
    midpointDefect (4/5) = 8/21 := by
  norm_num [midpointDefect]

end GNO
