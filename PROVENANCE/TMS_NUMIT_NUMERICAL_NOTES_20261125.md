# NuMIT sensitivity numerical checks

The first run stopped before treatment-effect calculations during independent
VAR calibration. The scipy root/Lyapunov solution matched its requested TDMI,
but evaluating a raw, poorly conditioned covariance with the imported explicit
cofactor determinant incurred cancellation (a reproduced draw differed by
1.55e-7 nats; the independent LAPACK log determinant differed by 1.49e-11).
Scaling this covariance to unit diagonal reduced the discrepancy to 5.44e-9.

The new calibration wrapper now checks the independently computed log determinant
and converts the covariance to a correlation matrix before analytic PhiID,
as the reference-null helper already does. This is an information-invariant
numerical conditioning step. No observed outcome, null prior, input data,
production helper, participant, or treatment comparison was changed. The
completed reference grids are reused, and calibration is rerun before effects.

A second synthetic covariance still exceeded 1e-6 after diagonal scaling.
Calibration observations therefore now use independent LAPACK log determinants
for all Gaussian information quantities. No synthetic draw is excluded for
cofactor inaccuracy; the number and largest discrepancy are saved. This makes
the population calibration a check on the reference helper's approximate
conditional distribution, including its numerical error, rather than comparing
two copies of the same cofactor calculation. The same PhiID lattice and MMI
definition are retained. Production TMS measures already matched the original
MATLAB results within the stated tolerance and are not changed.
