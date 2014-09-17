
# =========================================================================
#       Copyright (c) 2008 Stephen Lihn (stevelihn@gmail.com)
# =========================================================================
# Both this software and the term "Skew Lognormal Cascade Distribution" are
# registered with US Copyright Office. You can use this software freely 
# for personal and educational purpose as long as you retain this copyright 
# section and the DISCLAIMER OF WARRANTY in your application.
#
# Do not use this software and the term "Skew Lognormal Cascade Distribution"
# for commercial purpose without permission. 
# ==========================================================================
#                     DISCLAIMER OF WARRANTY. 
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTY OF ANY 
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, WITHOUT LIMITATION, 
# WARRANTIES THAT THE COVERED CODE IS FREE OF DEFECTS, MERCHANTABLE,
# FIT FOR A PARTICULAR PURPOSE OR NON-INFRINGING. THE ENTIRE RISK AS TO THE
# QUALITY AND PERFORMANCE OF THIS SOFTWARE IS WITH YOU. SHOULD THIS SOFTWARE
# PROVE DEFECTIVE IN ANY RESPECT, YOU ASSUME THE COST OF ANY NECESSARY 
# SERVICING, REPAIR OR CORRECTION. THIS DISCLAIMER OF WARRANTY CONSTITUTES 
# AN ESSENTIAL PART OF THE PERMISSION TO USE THIS SOFTWARE. NO USE 
# IS AUTHORIZED HEREUNDER EXCEPT UNDER THIS DISCLAIMER.
# ==========================================================================

# testing the moments of higher order lognormal cascade distribution
# whether the numbers match between direct integrals and analytical formula

more off;

beta = -0.15
gr   = 0.08
phi  = 1.6;

eta2  = 0.05
beta2 = 0.05
phi2  = 0.5

r2  = beta2*beta2*eta2*eta2 + 1.0;
eta = phi2*e^(eta2^2)*sqrt(r2)

# n=10, nsigma=20 is good for DOW parameters, eta=0.75
#       slog_dist_so_calc_stats ( -0.1, 0.08, 1.6, 0.1, 0.0, 0.7, 10, 20)

n=15;
nsigma1=30;
nsigma2=30;

  y1 = zeros(1,1);
  y2 = zeros(1,1);
  y3 = zeros(1,1);
  xrng  = floor(80*phi);
  xpart = xrng*20
  x = (-xpart:xpart)/xpart*xrng;
  dx = 1/xpart*xrng;

  y1 = slog_dist_simpson (x, eta, beta, gr, phi, 10, 10);
  len = length(y1)
  [ x(1), y1(1), x(len), y1(len) ]
  total_prob1 = sum(y1)*dx

  # y2 = slog_dist_so_simpson(x, beta, gr, phi, eta2, beta2, phi2, n, nsigma1 );

  xi = xpart+1;
  y2(xi) = slog_dist_so_simpson(x(xi), beta, gr, phi, eta2, beta2, phi2, n, nsigma1 );
  y2(2*xpart+1) = 0;

  st_theory = slog_dist_so_calc_stats ( beta, gr, phi, eta2, beta2, phi2, n, nsigma2 )
  
for i = 1 : xpart
  xi = xpart+1-i;
  y2(xi) = slog_dist_so_simpson(x(xi), beta, gr, phi, eta2, beta2, phi2, n, nsigma1 );
  xi = xpart+1+i;
  y2(xi) = slog_dist_so_simpson(x(xi), beta, gr, phi, eta2, beta2, phi2, n, nsigma1 );

if ( rem(i,10) == 0 || i == xpart )  

    p = [ i, xi, x(xi), log(y2(xi)), log(x(xi)^4*y2(xi)) ]
    total_prob2 = sum(y2)*dx

    u1 = sum(x.*y2)*dx;
    u2 = sum(x.^2.*y2)*dx;
    u3 = sum(x.^3.*y2)*dx;
    u4 = sum(x.^4.*y2)*dx;
  
    k1 = u1;
    f2_st(1) = u1;
    k2 = u2-u1*u1;
    f2_st(2) = sqrt(k2); # std = sqrt(variance)
    k3 = u3-3.0*k2*k1-k1*k1*k1;
    f2_st(3) = k3/sqrt(k2*k2*k2);
    k4 = u4-4.0*k3*k1-3.0*k2*k2-6.0*k2*k1*k1-k1*k1*k1*k1;
    f2_st(4) = k4/k2/k2;

    st_theory
    f2_st
endif

endfor
    