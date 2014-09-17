
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

# testing the integral of skew lognormal distribution

more off;

eta  = 0.7
beta = -0.1
gr   = 0.1
phi  = 1.6;

eta2  = 0.4
beta2 = 0.1
phi2  = eta/e^(eta2^2)
# phi2 can not be zero in SO, it should revert to first order 

y1 = zeros(1,1);
y2 = zeros(1,1);
y3 = zeros(1,1);
xrng  = 30;
xpart = xrng*5;
x = (-xpart:xpart)/xpart*xrng;

  y1 = slog_dist_simpson (x, eta, beta, gr, phi, 10, 10);
  y2 = slog_dist_so_simpson(x, beta, gr, phi, eta2, beta2, phi2, 20, 40);

  len = length(y1)
  [ x(1), y1(1) ]
  [ x(len), y1(len) ]
  
figure(1);
subplot(2,2,1);
plot(x, log(y1), x, log(y2));

# calculate the error between simple and simpson methods
subplot(2,2,2);
d2 = (log(y2)-log(y1))./log(y1);
plot(x, d2.*10^12);


# calculate total probability, should sum to one
total_prob1 = sum(y1)*1/xpart*xrng
total_prob2 = sum(y2)*1/xpart*xrng
prob_diff   = total_prob2 - total_prob1

# check convergence


xt = [ 0, 10, 20, 40, 80, 160, 320, 640, 1280 ];
eta2 =0.2
e.^xt .* slog_dist_simpson( xt, eta2, beta2, -beta2*eta2^2, phi2, 10, 10)
