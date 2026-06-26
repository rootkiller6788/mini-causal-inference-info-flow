# Pearl Structural Causal Models - Theory

## Structural Causal Model
A SCM M = (U, V, F, P(u)) consists of:
- Exogenous variables U (background, unobserved)
- Endogenous variables V (observed)
- Structural equations: v_i = f_i(PA_i, u_i)
- Distribution P(u) over exogenous variables

## Causal Graph
The directed graph G where nodes are V and edges PA_i -> V_i for each structural equation. Must be acyclic (DAG).

## The Ladder of Causation
1. Association: P(y|x) - seeing, observational
2. Intervention: P(y|do(x)) - doing, experimental
3. Counterfactuals: Y_x(u) - imagining, retrospective

## d-Separation (Pearl 1988)
A path p is blocked by Z if:
- p contains chain i->m->j or fork i<-m->j and m in Z, OR
- p contains collider i->m<-j and neither m nor its descendants in Z

## Back-Door Criterion
Z satisfies back-door relative to (X,Y) if:
1. No node in Z is a descendant of X
2. Z blocks all back-door paths from X to Y

Then: P(y|do(x)) = sum_z P(y|x,z)P(z)

## Front-Door Criterion
When back-door is blocked by unobserved confounders:
P(y|do(x)) = sum_m P(m|x) sum_x' P(y|x',m)P(x')
