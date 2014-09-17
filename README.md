Skew Lognormal Cascade Distribution
===================================

Implementation of the skew lognormal cascade distribution

# Status

Currently v1.2 of octave program is checked in under the /octave directory.
The plan is to rewrite it into either R package and/or C package (e.g. for possible inclusion into GNU GSL)

http://www.skew-lognormal-cascade-distribution.org/impl/

# Paper

The numerical method was published in SSRN: http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1273087

# Running Tests in Octave

* Run `slog_dist_test1.m` to watch the plot of the pdf. You can edit the parameters (ita, beta, gr) in the file to see how the distribution evolves. This test takes 10-20 seconds to run on a moderate PC.

* `slog_solve_zp_test.m` has more info about the inner structure of the integrand that leads to the pdf. 
