# Causal Effect Estimation with DAGs — Theory

## Structural Causal Models
An SCM M = (U, V, F) consists of exogenous U, endogenous V, and structural equations F.
Each v_i = f_i(pa_i, u_i) defines the causal mechanism.

## The do-operator
P(y | do(x)) is the post-intervention distribution after setting X=x via external intervention.
do(x) corresponds to removing the equation for X and setting it to x.

## Back-Door Criterion
A set Z satisfies the back-door criterion for (X, Y) iff:
1. No node in Z is a descendant of X
2. Z blocks every path between X and Y that starts with an arrow into X
Then: P(y | do(x)) = sum_z P(y | x, z) P(z)

## Front-Door Criterion
A mediator M satisfies the front-door criterion iff:
1. M intercepts all directed paths from X to Y
2. No unblocked back-door path X->...->M
3. All back-door paths M->Y are blocked by X
Then: P(y | do(x)) = sum_m P(m | x) sum_x' P(y | x', m) P(x')

## Estimation Methods
- IPW: reweight by 1/P(T|X)
- G-Computation: standardize by outcome model
- Doubly Robust: combine IPW + outcome model, consistent if either is correct

References: Pearl (2009), Hernan & Robins (2020), Rubin (1974)
