import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Tactic

namespace GNO

/-- The metric core of the paper's two-point minimax lower bound. -/
theorem two_point_lower_bound {X : Type*} [PseudoMetricSpace X]
    (x y estimate : X) (r : ℝ) (hsep : 2*r ≤ dist x y) :
    r ≤ max (dist x estimate) (dist y estimate) := by
  have ht := dist_triangle x estimate y
  rw [dist_comm estimate y] at ht
  have h1 := le_max_left (dist x estimate) (dist y estimate)
  have h2 := le_max_right (dist x estimate) (dist y estimate)
  linarith

/-- Pure factorization. For a stochastic-kernel application, Y can be a
type of probability measures; measurability is a separate obligation. -/
theorem factors_iff_fiber_constant {X C Y : Type*}
    (observe : X → C) (K : X → Y) (hsurj : Function.Surjective observe) :
    (∃ coarse : C → Y, ∀ x, coarse (observe x) = K x) ↔
    (∀ x y, observe x = observe y → K x = K y) := by
  constructor
  · rintro ⟨coarse, hc⟩ x y hxy
    rw [← hc x, ← hc y, hxy]
  · intro hfiber
    classical
    let representative : C → X := fun c => Classical.choose (hsurj c)
    have hrep : ∀ c, observe (representative c) = c :=
      fun c => Classical.choose_spec (hsurj c)
    refine ⟨fun c => K (representative c), ?_⟩
    intro x
    exact hfiber (representative (observe x)) x (hrep (observe x))

end GNO
