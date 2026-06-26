# Granger Causality Theory

## Definition (Granger 1969)
A time series X Granger-causes Y if past values of X contain
information that helps predict Y beyond what is contained in
past values of Y alone.

## VAR Formulation
Restricted:  Y_t = c + sum(alpha_i * Y_{t-i}) + eps_t
Unrestricted: Y_t = c + sum(alpha_i * Y_{t-i}) + sum(beta_i * X_{t-i}) + eta_t

## F-test
F = ((RSS_r - RSS_ur) / p) / (RSS_ur / (n - 2p - 1))

H0: beta_1 = ... = beta_p = 0 (X does NOT Granger-cause Y)
Reject H0 if F > F_critical or p < alpha

## References
Granger, C.W.J. (1969). Econometrica, 37(3):424-438.
Geweke, J. (1982). JASA, 77(378):304-313.
