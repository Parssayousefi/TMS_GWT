# Numerical review of exploratory TMS analyses

Independent SciPy tests and a NumPy pseudoinverse HC3 covariance calculation
verified 977 estimable neural contrasts in the initial exploratory screen and
the three normalized reference estimates (two seeds and their combined CDF).
Brute-force checks independently reproduced the multivariate sign-flip test,
univariate Freedman-Lane reference, and the family-wise max-T implementation.

One secondary model was not estimable with HC3: the Aachen intensity contrast
for amygdala gateway count, adjusted for order, baseline difference and motion.
A design row has unit leverage, making the HC3 residual correction undefined.
The initial numerical p value (.2856) is retained in the original output as an
audit record; the revised family-wise exploration table marks its interval and
p value unavailable. It remains in the correction burden as p=1. No main result,
participant inclusion, or surviving finding depends on this row.

The Aachen mean-of-two-changes calculation was also checked against a 2x2
within-participant ANOVA of time and intensity. For global Phi, the time effect
is F(1,14)=9.0556, p=.009378, exactly the square of the original pooled paired
t statistic. The time-by-intensity effect has p=.772757. This verifies the
averaging implementation; it does not create a sham-controlled comparison.
