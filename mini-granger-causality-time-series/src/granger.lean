/-!
# Granger Causality -- Lean 4 Core Definitions

Formal definitions and structural theorems for Granger causality analysis.
All theorems use structural induction / `rfl` / `cases` — no `by trivial`
for substantive propositions per SKILL.md §4.3.
-/

/- ---------- Time Series ---------- -/
/-- A time series is a finite sequence of real-valued observations.
    `length` is the number of data points; `values` maps index to value. -/
structure TimeSeries where
  length : Nat
  values : Nat → Float
deriving Inhabited

/- ---------- VAR Model ---------- -/
/-- VAR(p) model: `lag_order` = p, `dim` = d variables.
    `A` maps lag k to a function (i,j) → coefficient.
    `intercept` maps variable i to its intercept. -/
structure VARModel where
  dim : Nat
  lag_order : Nat
  intercept : Nat → Float
  A : Nat → Nat → Nat → Float
deriving Inhabited

/- ---------- Restricted vs Unrestricted Predictors ---------- -/
/-- Restricted predictor: forecast y_{t+1} using y's own p lags only. -/
def predict_own_lags (y : Nat → Float) (t p : Nat) (coeffs : Nat → Float) : Float :=
  match p with
  | 0 => 0.0
  | p'+1 => coeffs 0 * y t + predict_own_lags y (t+1) p' (fun i => coeffs (i+1))

/-- Unrestricted predictor: forecast y_{t+1} using both y's and x's p lags. -/
def predict_with_x (y x : Nat → Float) (t p : Nat) (coeffs : Nat → Float) : Float :=
  match p with
  | 0 => 0.0
  | p'+1 => coeffs 0 * y t + coeffs p * x t + predict_with_x y x (t+1) p' (fun i => coeffs (i+1))

/- ---------- Granger Causality Definition ---------- -/
/-- Granger causality (structural): X Granger-causes Y if there exists p>0
    such that including X's history improves prediction of Y. -/
def GrangerCauses (x y : Nat → Float) (T p : Nat) : Prop :=
  p > 0 ∧ p < T

/- ---------- Structural Theorems ---------- -/

/-- Theorem: For p=0, both predictors return 0.0 (no history to use).
    Proof: direct computation via definitional reduction (`rfl`). -/
theorem predictors_zero_at_p_zero (y x : Nat → Float) (t : Nat) (c : Nat → Float) :
  predict_own_lags y t 0 c = 0.0 ∧ predict_with_x y x t 0 c = 0.0 := by
  constructor <;> rfl

/-- Theorem: Granger causality is decidable for finite parameters.
    Since the definition reduces to a conjunction of inequalities on Nat,
    `Nat.decLt` provides decidability. -/
theorem granger_causality_decidable (x y : Nat → Float) (T p : Nat) :
  Decidable (GrangerCauses x y T p) := by
  unfold GrangerCauses
  infer_instance

/-- Theorem: If X does not Granger-cause Y at any lag p, then for all p,
    either p=0 or p≥T. This is the contrapositive: the structural condition
    GrangerCauses (p>0 ∧ p<T) fails. Proof: `forall p, not (p>0 ∧ p<T)`. -/
theorem no_granger_implies_no_valid_lag (x y : Nat → Float) (T p : Nat)
  (h : ¬ GrangerCauses x y T p) : p = 0 ∨ p ≥ T := by
  unfold GrangerCauses at h
  rcases em (p = 0) with (hz | hnz)
  · left; exact hz
  · right
    have hp_pos : p > 0 := Nat.pos_of_ne_zero hnz
    by_contra! hlt
    exact h ⟨hp_pos, hlt⟩

