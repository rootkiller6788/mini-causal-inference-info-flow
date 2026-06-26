# Transfer Entropy

## Definition (Schreiber 2000)
Transfer entropy from X to Y quantifies the reduction in
uncertainty about future values of Y given past values of
X, beyond the information already contained in past values of Y.

TE_{X->Y} = I(Y_{t+1}; X_t | Y_t)
          = H(Y_{t+1}|Y_t) - H(Y_{t+1}|Y_t, X_t)

## Relationship to Granger Causality
For Gaussian variables, transfer entropy is equivalent to
Granger causality (Barnett et al. 2009).

For non-Gaussian/nonlinear systems, transfer entropy captures
information flow that linear Granger causality would miss.

## References
Schreiber, T. (2000). Physical Review Letters, 85(2):461.
Barnett, L. et al. (2009). Physical Review Letters, 103(23):238701.
