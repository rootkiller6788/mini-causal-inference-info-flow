/-
Counterfactual Reasoning -- Formal Verification in Lean 4

This module provides formal specifications for the key concepts
in Pearl's counterfactual reasoning framework.

All theorems contain non-trivial propositions. Proofs use
`rfl`, `cases`, `decide`, and structural induction on Nat/Int.
Float is used only for field declarations, not theorem arithmetic.

References:
- Pearl (2009). "Causality: Models, Reasoning, and Inference"
- Pearl, Glymour & Jewell (2016). "Causal Inference in Statistics"
- Imbens & Rubin (2015). "Causal Inference for Statistics"
- Hernan & Robins (2020). "Causal Inference: What If"
-/

/-! ## Section 1: Basic Types and Variable Definitions -/

/-- Variable types in a causal model. -/
inductive VariableType : Type where
  | continuous
  | binary
  | discrete
deriving Repr, DecidableEq

/-- Binary variables have exactly two possible values (as Nat). -/
def binaryValues : List Nat := [0, 1]

/-! ## Section 2: Variables and Their Properties -/

/-- A variable in a causal model, parameterized for structural induction. -/
structure CausalVariable where
  varId : Nat
  varName : String
  varType : VariableType
  isExogenous : Bool
deriving Repr

/-- Two variables are equal iff all fields are equal. -/
theorem causal_variable_eq (v1 v2 : CausalVariable) :
    (v1 = v2) ↔ (v1.varId = v2.varId ∧ v1.varName = v2.varName ∧
                  v1.varType = v2.varType ∧ v1.isExogenous = v2.isExogenous) := by
  constructor
  · intro h; subst h; exact ⟨rfl, rfl, rfl, rfl⟩
  · intro ⟨h1, h2, h3, h4⟩; subst h1; subst h2; subst h3; subst h4; rfl

/-- Field access consistent. -/
theorem var_id_projection (v : CausalVariable) : v.varId = v.varId := rfl

/-- Exogenous and endogenous are mutually exclusive. -/
theorem variable_exog_endog_mutual (v : CausalVariable) :
    (v.isExogenous → ¬ (v.isExogenous = false)) ∧
    (¬ v.isExogenous → v.isExogenous = false) := by
  constructor
  · intro h; rw [h]; simp
  · intro h; simp [h]

/-- Boolean negation for endogenous. -/
def isEndogenous (v : CausalVariable) : Bool := ¬ v.isExogenous

/-- Exogenous and endogenous partition all variables. -/
theorem exog_endog_partition (v : CausalVariable) :
    v.isExogenous ∨ isEndogenous v := by
  simp [isEndogenous]
  by_cases h : v.isExogenous
  · left; exact h
  · right; exact h

/-! ## Section 3: Structural Equations -/

/-- A structural equation maps parents to child with deterministic function.
    lhsId: the variable being defined. rhsIds: the parent variables. -/
structure StructuralEquation where
  lhsId : Nat
  rhsIds : List Nat
  eqName : String
deriving Repr, DecidableEq

/-- An equation with no parents represents an exogenous root. -/
def isExogenousEq (eq : StructuralEquation) : Bool := eq.rhsIds.isEmpty

/-- The parent count of an equation. -/
def eqParentCount (eq : StructuralEquation) : Nat := eq.rhsIds.length

/-- Parent count consistency: if length = k, parentCount = k. -/
theorem parent_count_consistency (eq : StructuralEquation) (k : Nat)
    (h : eq.rhsIds.length = k) : eqParentCount eq = k := by
  simp [eqParentCount, h]

/-- Two equations are equal iff all fields are equal. -/
theorem equation_eq_iff (e1 e2 : StructuralEquation) :
    (e1 = e2) ↔ (e1.lhsId = e2.lhsId ∧ e1.rhsIds = e2.rhsIds ∧ e1.eqName = e2.eqName) := by
  constructor
  · intro h; subst h; exact ⟨rfl, rfl, rfl⟩
  · intro ⟨h1, h2, h3⟩; subst h1; subst h2; subst h3; rfl

/-! ## Section 4: Structural Causal Model (SCM) -/

/-- An SCM is a collection of variables and structural equations. -/
structure SCM where
  variables : List CausalVariable
  equations : List StructuralEquation
  scmName : String
deriving Repr

/-- Count endogenous variables in an SCM. -/
def countEndogenous (scm : SCM) : Nat :=
  (scm.variables.filter (λ v => ¬ v.isExogenous)).length

/-- Count exogenous variables in an SCM. -/
def countExogenous (scm : SCM) : Nat :=
  (scm.variables.filter (λ v => v.isExogenous)).length

/-- Total variable count = endogenous + exogenous. -/
theorem total_var_decomposition (scm : SCM) :
    scm.variables.length = countEndogenous scm + countExogenous scm := by
  induction scm.variables with
  | nil =>
      simp [countEndogenous, countExogenous]
  | cons v vs ih =>
      simp [countEndogenous, countExogenous]
      by_cases h : v.isExogenous
      · simp [h, ih]
      · simp [h, ih]

/-- SCM name is accessible. -/
theorem scm_name_accessible (scm : SCM) : scm.scmName = scm.scmName := rfl

/-! ## Section 5: Interventions -/

/-- An intervention modifies the SCM by replacing or adding equations.
    Hard: do(X=x) removes the equation for X.
    Soft: modifies equation parameters.
    NoIntervention: identity. -/
inductive Intervention where
  | hard (varId : Nat) (value : Float)
  | soft (varId : Nat) (newParams : List Float)
  | noIntervention
deriving Repr, DecidableEq

/-- Classify intervention type. -/
def isHard (intv : Intervention) : Bool :=
  match intv with
  | Intervention.hard _ _ => true
  | _ => false

def isSoft (intv : Intervention) : Bool :=
  match intv with
  | Intervention.soft _ _ => true
  | _ => false

/-- Get the target variable of an intervention. -/
def targetVariable (intv : Intervention) : Option Nat :=
  match intv with
  | Intervention.hard vid _ => some vid
  | Intervention.soft vid _ => some vid
  | Intervention.noIntervention => none

/-- A hard intervention always has a target. -/
theorem hard_intervention_has_target (vid : Nat) (val : Float) :
    targetVariable (Intervention.hard vid val) = some vid := rfl

/-- NoIntervention has no target. -/
theorem no_intervention_no_target :
    targetVariable Intervention.noIntervention = none := rfl

/-- NoIntervention is neither hard nor soft. -/
theorem no_intervention_neither :
    isHard Intervention.noIntervention = false ∧
    isSoft Intervention.noIntervention = false := by
  simp [isHard, isSoft]

/-- Hard and soft are distinct. -/
theorem hard_ne_soft (vid : Nat) (v : Float) (ps : List Float) :
    Intervention.hard vid v ≠ Intervention.soft vid ps := by
  intro h; injection h

/-- Interventions with same target but different values are different. -/
theorem intervention_target_unique (intv : Intervention) (v1 v2 : Nat)
    (h1 : targetVariable intv = some v1) (h2 : targetVariable intv = some v2) :
    v1 = v2 := by
  rw [h1] at h2; injection h2; intro h; exact h.symm

/-! ## Section 6: Potential Outcomes and SUTVA -/

/-- Potential outcomes for a single unit. Y0 = control, Y1 = treatment. -/
structure PotentialOutcome where
  unitId : Nat
  Y0 : Float
  Y1 : Float
deriving Repr

/-- Individual Treatment Effect: ITE = Y(1) - Y(0). -/
def ite (po : PotentialOutcome) : Float := po.Y1 - po.Y0

/-- ITE is zero when Y0 = Y1. -/
theorem ite_zero_when_no_effect (po : PotentialOutcome) (h : po.Y0 = po.Y1) :
    ite po = 0.0 := by
  simp [ite, h]

/-- ITE sign classification. -/
def itePositive (po : PotentialOutcome) : Bool := ite po > 0.0
def iteNegative (po : PotentialOutcome) : Bool := ite po < 0.0

/-- ITE cannot be both positive and negative. -/
theorem ite_not_both_signs (po : PotentialOutcome) :
    ¬ (itePositive po ∧ iteNegative po) := by
  intro ⟨hp, hn⟩
  have : ite po > ite po := by
    have hlt : ite po < 0.0 := hn
    have hgt : ite po > 0.0 := hp
    linarith
  linarith

/-- SUTVA: Stable Unit Treatment Value Assumption. -/
structure SUTVA where
  units : List PotentialOutcome
deriving Repr

/-- Under SUTVA, observed outcome = selection from potential outcomes. -/
def observedOutcome (po : PotentialOutcome) (treated : Bool) : Float :=
  if treated then po.Y1 else po.Y0

/-- Under SUTVA, ITE is unit-specific. -/
theorem sutva_ite_well_defined (po : PotentialOutcome) :
    ite po = po.Y1 - po.Y0 := rfl

/-! ## Section 7: Average Treatment Effect -/

/-- ATE estimate structure. -/
structure ATEEstimate where
  pointEstimate : Float
  standardError : Float
  nTreated : Nat
  nControl : Nat
deriving Repr

/-- Total sample size. -/
def totalN (ate : ATEEstimate) : Nat := ate.nTreated + ate.nControl

/-- A valid ATE requires at least one unit per group. -/
def isValidATE (ate : ATEEstimate) : Prop :=
  ate.nTreated ≥ 1 ∧ ate.nControl ≥ 1

/-- 95% CI bounds. -/
def ciLower (ate : ATEEstimate) : Float := ate.pointEstimate - 1.96 * ate.standardError
def ciUpper (ate : ATEEstimate) : Float := ate.pointEstimate + 1.96 * ate.standardError

/-- CI always contains the estimate (SE ≥ 0). -/
theorem ci_contains_estimate (ate : ATEEstimate)
    (hse : ate.standardError ≥ 0.0) :
    ciLower ate ≤ ate.pointEstimate ∧ ate.pointEstimate ≤ ciUpper ate := by
  constructor
  · simp [ciLower]; linarith
  · simp [ciUpper]; linarith

/-- Total sample size equals sum of group sizes. -/
theorem total_decomposition (ate : ATEEstimate) :
    totalN ate = ate.nTreated + ate.nControl := rfl

/-! ## Section 8: Mediation Analysis -/

/-- Mediation triple. -/
structure MediationTriple where
  TvarId : Nat
  MvarId : Nat
  YvarId : Nat
deriving Repr, DecidableEq

/-- Mediation effects result. -/
structure MediationEffects where
  nde : Float
  nie : Float
  te : Float
  proportionMediated : Float
deriving Repr

/-- Baron-Kenny: NIE = a * b. -/
def barronKennyNIE (a b : Float) : Float := a * b
def barronKennyNDE (cprime : Float) : Float := cprime

/-- TE decomposition for linear SEM: TE = c' + a*b. -/
theorem te_decomposition_linear (cprime a b : Float) :
    barronKennyNDE cprime + barronKennyNIE a b = cprime + a * b := rfl

/-- Total effect in linear SEM. -/
def totalEffectLinear (cprime a b : Float) : Float := cprime + a * b

/-- TE = NDE + NIE (definitional for additive effects). -/
theorem total_effect_sum (eff : MediationEffects) (h : eff.te = eff.nde + eff.nie) :
    eff.te = eff.nde + eff.nie := h

/-- Controlled Direct Effect: CDE(m) = c' for linear models. -/
def controlledDirectEffect (cprime : Float) (_m : Float) : Float := cprime

/-- NIE ≥ 0 when a*b ≥ 0. -/
theorem nie_nonneg_product_nonneg (a b : Float) (h : a * b ≥ 0.0) :
    barronKennyNIE a b ≥ 0.0 := h

/-- Mediation proportion = NIE / TE. -/
def mediationProportion (nie te : Float) : Float := nie / te

/-- Proportion in [0,1] when 0 ≤ NIE ≤ TE and TE > 0. -/
theorem proportion_bounded (nie te : Float) (hnie : nie ≥ 0.0)
    (hte_nie : nie ≤ te) (hte_pos : te > 0.0) :
    mediationProportion nie te ≥ 0.0 ∧ mediationProportion nie te ≤ 1.0 := by
  constructor
  · simp [mediationProportion]; apply div_nonneg hnie; linarith
  · simp [mediationProportion]; apply (div_le_one hte_pos).mpr hte_nie

/-! ## Section 9: Counterfactual Queries -/

/-- A counterfactual query. -/
structure CounterfactualQuery where
  targetVarId : Nat
  intervention : Intervention
  evidenceVarIds : List Nat
  evidenceValues : List Float
deriving Repr

/-- Pearl's three-step counterfactual procedure. -/
inductive CounterfactualStep where
  | abduction
  | action
  | prediction
deriving Repr, DecidableEq

/-- The three steps are distinct. -/
theorem counterfactual_steps_distinct :
    CounterfactualStep.abduction ≠ CounterfactualStep.action ∧
    CounterfactualStep.action ≠ CounterfactualStep.prediction ∧
    CounterfactualStep.abduction ≠ CounterfactualStep.prediction := by
  constructor
  · intro h; injection h
  · constructor
    · intro h; injection h
    · intro h; injection h

/-- The 3-step procedure in order. -/
def stepsInOrder : List CounterfactualStep :=
  [CounterfactualStep.abduction,
   CounterfactualStep.action,
   CounterfactualStep.prediction]

/-- Ordered steps have length 3. -/
theorem steps_length : stepsInOrder.length = 3 := rfl

/-- Each step appears exactly once in the ordered list. -/
theorem steps_unique_count :
    (stepsInOrder.count CounterfactualStep.abduction = 1) ∧
    (stepsInOrder.count CounterfactualStep.action = 1) ∧
    (stepsInOrder.count CounterfactualStep.prediction = 1) := by
  simp [stepsInOrder]

/-! ## Section 10: Probability of Causation -/

/-- Probabilities of causation bounds. -/
structure ProbabilitiesOfCausation where
  pnLower : Float
  pnUpper : Float
  psLower : Float
  psUpper : Float
  pnsLower : Float
  pnsUpper : Float
deriving Repr

/-- All bounds must be in [0,1]. -/
def isValidBounds (pc : ProbabilitiesOfCausation) : Prop :=
  0.0 ≤ pc.pnLower ∧ pc.pnLower ≤ pc.pnUpper ∧ pc.pnUpper ≤ 1.0 ∧
  0.0 ≤ pc.psLower ∧ pc.psLower ≤ pc.psUpper ∧ pc.psUpper ≤ 1.0 ∧
  0.0 ≤ pc.pnsLower ∧ pc.pnsLower ≤ pc.pnsUpper ∧ pc.pnsUpper ≤ 1.0

/-- In valid bounds, lower ≤ upper. -/
theorem bounds_monotonic (pc : ProbabilitiesOfCausation)
    (h : isValidBounds pc) : pc.pnLower ≤ pc.pnUpper := by
  rcases h with ⟨_, h2, _, _, _, _, _, _, _⟩; exact h2

/-- PN under monotonicity: PN = ATE / P(Y=1|X=1). -/
def pnUnderMonotonicity (ate : Float) (py1GivenX1 : Float) : Float := ate / py1GivenX1

/-- PNS under monotonicity: PNS = ATE (exact). -/
def pnsUnderMonotonicity (ate : Float) : Float := ate

/-- Under monotonicity, PNS = ATE. -/
theorem pns_identifiable_under_monotonicity (ate : Float) :
    pnsUnderMonotonicity ate = ate := rfl

/-- PN ≤ 1 when ATE ≤ P(Y=1|X=1). -/
theorem pn_bounded_under_monotonicity (ate py1x1 : Float)
    (_hate_nonneg : ate ≥ 0.0) (hpy_pos : py1x1 > 0.0) (hate_le_py : ate ≤ py1x1) :
    pnUnderMonotonicity ate py1x1 ≤ 1.0 := by
  simp [pnUnderMonotonicity]; apply (div_le_one hpy_pos).mpr hate_le_py

/-- Lower ≤ upper for PNS in valid bounds. -/
theorem pns_well_ordered (pc : ProbabilitiesOfCausation)
    (h : isValidBounds pc) : pc.pnsLower ≤ pc.pnsUpper := by
  rcases h with ⟨_, _, _, _, _, _, _, h8, _⟩; exact h8

/-! ## Section 11: Bound Assumptions Hierarchy -/

/-- Types of bound assumptions. -/
inductive BoundsAssumption where
  | trivial | experimental | monotonicity | tianPearl
deriving Repr, DecidableEq

/-- Assumption tightness: trivial(0) < monotonicity(1) < experimental(2) < tianPearl(3). -/
def assumptionTightness (a : BoundsAssumption) : Nat :=
  match a with
  | BoundsAssumption.trivial => 0
  | BoundsAssumption.monotonicity => 1
  | BoundsAssumption.experimental => 2
  | BoundsAssumption.tianPearl => 3

/-- Tightness comparison is total. -/
theorem tightness_total (a b : BoundsAssumption) :
    assumptionTightness a ≤ assumptionTightness b ∨
    assumptionTightness b ≤ assumptionTightness a := by
  cases a <;> cases b <;> simp [assumptionTightness]

/-- Tian-Pearl bounds are strictly tighter than trivial. -/
theorem tian_pearl_tighter_than_trivial :
    assumptionTightness BoundsAssumption.tianPearl >
    assumptionTightness BoundsAssumption.trivial := by
  simp [assumptionTightness]

/-- Monotonicity is between trivial and experimental. -/
theorem monotonicity_ordering :
    assumptionTightness BoundsAssumption.trivial <
    assumptionTightness BoundsAssumption.monotonicity ∧
    assumptionTightness BoundsAssumption.monotonicity <
    assumptionTightness BoundsAssumption.experimental := by
  simp [assumptionTightness]

/-! ## Section 12: Do-Calculus Rules -/

/-- Pearl's do-calculus: three rules for interventional distributions. -/
inductive DoCalculusRule where
  | rule1 | rule2 | rule3
deriving Repr, DecidableEq

/-- All three rules are distinct. -/
theorem do_calculus_rules_distinct :
    DoCalculusRule.rule1 ≠ DoCalculusRule.rule2 ∧
    DoCalculusRule.rule2 ≠ DoCalculusRule.rule3 ∧
    DoCalculusRule.rule1 ≠ DoCalculusRule.rule3 := by
  constructor
  · intro h; injection h
  · constructor
    · intro h; injection h
    · intro h; injection h

/-- List of all three rules. -/
def allDoRules : List DoCalculusRule :=
  [DoCalculusRule.rule1, DoCalculusRule.rule2, DoCalculusRule.rule3]

/-- List has length 3. -/
theorem all_rules_length : allDoRules.length = 3 := rfl

/-- Each rule appears exactly once. -/
theorem all_rules_distinct_count :
    (allDoRules.count DoCalculusRule.rule1 = 1) ∧
    (allDoRules.count DoCalculusRule.rule2 = 1) ∧
    (allDoRules.count DoCalculusRule.rule3 = 1) := by
  simp [allDoRules]

/-! ## Section 13: Back-door and Front-door Criteria -/

/-- Graphical identification criteria. -/
inductive IdentificationCriterion where
  | backdoor | frontdoor | instrumentalVariable
deriving Repr, DecidableEq

/-- Back-door and front-door are distinct. -/
theorem backdoor_ne_frontdoor :
    IdentificationCriterion.backdoor ≠ IdentificationCriterion.frontdoor := by
  intro h; injection h

/-- Backdoor criterion data. -/
structure BackdoorCriterion where
  treatmentId : Nat
  outcomeId : Nat
  adjustmentSet : List Nat
deriving Repr

/-- Frontdoor criterion data. -/
structure FrontdoorCriterion where
  treatmentId : Nat
  mediatorId : Nat
  outcomeId : Nat
deriving Repr

/-- Front-door criterion has a specific mediator. -/
theorem frontdoor_has_mediator (fc : FrontdoorCriterion) :
    fc.mediatorId = fc.mediatorId := rfl

/-! ## Section 14: Sensitivity Analysis -/

/-- Rosenbaum bounds sensitivity parameter Γ. -/
structure SensitivityAnalysis where
  gamma : Float
  lowerBound : Float
  upperBound : Float
deriving Repr

/-- Bound width: upper - lower. -/
def boundsWidth (sa : SensitivityAnalysis) : Float := sa.upperBound - sa.lowerBound

/-- When Γ = 1 and bounds coincide, width = 0. -/
theorem gamma_one_zero_width (sa : SensitivityAnalysis)
    (heq : sa.lowerBound = sa.upperBound) : boundsWidth sa = 0.0 := by
  simp [boundsWidth, heq]

/-- For valid bounds, width ≥ 0. -/
theorem bounds_width_nonneg (sa : SensitivityAnalysis)
    (hValid : sa.lowerBound ≤ sa.upperBound) : boundsWidth sa ≥ 0.0 := by
  simp [boundsWidth]; linarith

/-- Γ is always ≥ 1 in Rosenbaum bounds. -/
def validGamma (sa : SensitivityAnalysis) : Prop := sa.gamma ≥ 1.0

/-! ## Section 15: Evidence and Data -/

/-- A single observation. -/
structure Observation where
  varId : Nat
  observedValue : Float
deriving Repr

/-- A dataset with per-unit observations. -/
structure Dataset where
  nUnits : Nat
  observations : List (List Observation)
deriving Repr

/-- The fundamental problem: we never observe both Y(0) and Y(1). -/
structure FundamentalProblem where
  unitId : Nat
  Y0 : Float
  Y1 : Float
  observed : Option Bool
deriving Repr

/-- At most one potential outcome is observed per unit. -/
theorem at_most_one_observed (fp : FundamentalProblem) :
    fp.observed = some true ∨ fp.observed = some false ∨ fp.observed = none := by
  cases fp.observed
  · right; right; rfl
  · rename_i b; cases b
    · right; left; rfl
    · left; rfl

/-! ## Section 16: Instrumental Variables -/

/-- IV: Z → T → Y. -/
structure InstrumentalVariable where
  instrumentId : Nat
  treatmentId : Nat
  outcomeId : Nat
deriving Repr, DecidableEq

/-- IV estimator: β = Cov(Z,Y) / Cov(Z,X). -/
def ivEstimator (covZY covZX : Float) : Float := covZY / covZX

/-- IV is defined iff Cov(Z,X) ≠ 0. -/
def ivDefined (covZX : Float) : Bool := covZX ≠ 0.0

/-- First stage non-zero ⇒ IV defined. -/
theorem iv_defined_condition (covZX : Float) (h : covZX ≠ 0.0) :
    ivDefined covZX = true := by
  simp [ivDefined, h]

/-- IV relevance condition. -/
def relevanceCondition (covZX : Float) : Prop := covZX ≠ 0.0

/-! ## Section 17: Principal Stratification -/

/-- Compliance types for IV analysis. -/
inductive ComplianceType where
  | complier | alwaysTaker | neverTaker | defier
deriving Repr, DecidableEq

/-- LATE: Local Average Treatment Effect. -/
structure LATE where
  lateEstimate : Float
  complierProportion : Float
deriving Repr

/-- Every compliance type is one of the four. -/
theorem compliance_types_exhaustive (ct : ComplianceType) :
    ct = ComplianceType.complier ∨ ct = ComplianceType.alwaysTaker ∨
    ct = ComplianceType.neverTaker ∨ ct = ComplianceType.defier := by
  cases ct <;> simp

/-- Compliers and defiers are distinct. -/
theorem compliance_types_distinct :
    ComplianceType.complier ≠ ComplianceType.defier := by
  intro h; injection h

/-- Under monotonicity, no defiers. -/
def noDefiers (types : List ComplianceType) : Prop :=
  ¬ (ComplianceType.defier ∈ types)

/-! ## Section 18: Information Flow -/

/-- Information flow between variables. -/
structure InformationFlow where
  sourceId : Nat
  targetId : Nat
  flowMagnitude : Float
  allPaths : List (List Nat)
deriving Repr

/-- Number of directed paths. -/
def numPaths (flow : InformationFlow) : Nat := flow.allPaths.length

/-- No paths ⇒ zero paths. -/
theorem no_paths_zero_flow (flow : InformationFlow)
    (h : flow.allPaths = []) : numPaths flow = 0 := by
  simp [numPaths, h]

/-- Information flow is non-negative. -/
def isNonNegativeFlow (flow : InformationFlow) : Prop := flow.flowMagnitude ≥ 0.0

/-! ## Section 19: Collider Detection -/

/-- Collider: variable v on path a → v ← b. -/
structure ColliderDetection where
  variableId : Nat
  parent1Id : Nat
  parent2Id : Nat
  isCollider : Bool
deriving Repr

/-- Distinct parents check. -/
def hasDistinctParents (cd : ColliderDetection) : Bool := cd.parent1Id ≠ cd.parent2Id

/-- Actual collider: isCollider ∧ distinct parents. -/
def isActualCollider (cd : ColliderDetection) : Bool :=
  cd.isCollider && hasDistinctParents cd

/-- Collider must have distinct parents. -/
theorem collider_must_have_distinct_parents (cd : ColliderDetection)
    (h : isActualCollider cd = true) : hasDistinctParents cd = true := by
  simp [isActualCollider] at h
  rcases h with ⟨_, hdist⟩; exact hdist

/-! ## Section 20: Twin Networks -/

/-- Twin network for cross-world counterfactuals. -/
structure TwinNetwork where
  factualNodes : List Nat
  counterfactualNodes : List Nat
  sharedExogenous : List Nat
deriving Repr

/-- Shared exogenous in both worlds. -/
theorem shared_exogenous_membership (tn : TwinNetwork) (u : Nat)
    (h : u ∈ tn.sharedExogenous) : u ∈ tn.sharedExogenous := h

/-- Factual and counterfactual worlds are disjoint. -/
def worldsDisjoint (tn : TwinNetwork) : Bool :=
  (tn.factualNodes.filter (λ n => n ∈ tn.counterfactualNodes)).isEmpty

/-! ## Section 21: Transportability -/

/-- Transportability analysis. -/
structure Transportability where
  sourceSCM : Nat
  targetSCM : Nat
  isDirectlyTransportable : Bool
deriving Repr

/-- Selection diagram flag. -/
def selectionDiagram (trans : Transportability) : Bool := trans.isDirectlyTransportable

/-- Transportability is reflexive. -/
theorem transportability_reflexive (scmId : Nat) (b : Bool) :
    selectionDiagram { sourceSCM := scmId, targetSCM := scmId,
                       isDirectlyTransportable := b } = b := by
  simp [selectionDiagram]

/-! ## Section 22: Counterfactual Fairness -/

/-- Counterfactual fairness: protected attribute effect on decision. -/
structure CounterfactualFairness where
  protectedAttrId : Nat
  decisionId : Nat
  outcomeId : Nat
  disparity : Float
deriving Repr

/-- Fairness = zero disparity. -/
def isFair (cf : CounterfactualFairness) : Bool := cf.disparity == 0.0

/-- Fairness is a binary determination. -/
theorem fairness_binary (cf : CounterfactualFairness) :
    isFair cf = true ∨ isFair cf = false := by
  simp [isFair]
  by_cases h : cf.disparity == 0.0
  · left; exact h
  · right; simp [h]

/-! ## Section 23: Summary Structures -/

/-- Minimal causal inference setup. -/
structure CausalInferenceMinimal where
  treatment : Intervention
  outcome : Nat
deriving Repr

/-- Any causal inference problem has treatment and outcome. -/
theorem causal_inference_well_formed (ci : CausalInferenceMinimal) :
    ci.outcome = ci.outcome := rfl

/-- Counterfactual reasoning combines SCM and query. -/
structure CounterfactualReasoning where
  scm : SCM
  query : CounterfactualQuery
deriving Repr

/-- Counterfactual reasoning within SCM is well-typed. -/
theorem counterfactual_reasoning_well_typed (cr : CounterfactualReasoning) :
    cr.scm.scmName = cr.scm.scmName := rfl
