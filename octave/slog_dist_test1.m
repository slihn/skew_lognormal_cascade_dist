
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

ita  = 0.5
beta = 0.5
gr   = 0.5
phi  = 1;

y1 = zeros(1,1);
y2 = zeros(1,1);
y3 = zeros(1,1);
xrng  = 30;
xpart = xrng*30;
x = (-xpart:xpart)/xpart*xrng;

for xi = 1 : length(x)
  y1(xi) = slog_dist_simple (x(xi), ita, beta, gr, phi, 15, 5000);
  y2(xi) = slog_dist_simpson(x(xi), ita, beta, gr, phi, 10, 10);
  y3(xi) = slog_dist_simple (x(xi), ita, 0.0, 0.0, phi, 15, 5000);
endfor
  len = length(y1)
  [ x(1), y1(1) ]
  [ x(len), y1(len) ]
  
figure(1);
subplot(2,2,1);
plot(x, log(y1), x, log(y3), x, log(y2));

# calculate the error between simple and simpson methods
subplot(2,2,2);
d2 = (log(y2)-log(y1))./log(y1);
plot(x, d2.*10^12);

# diff between skew and symmetric curves
subplot(2,2,3);
d3 = (log(y3)-log(y1))./log(y1);
plot(x, d3.*10^6);

# test slope evolution at the tail
subplot(2,2,4);
xi2 = (length(y1)-xpart+40):length(y1);
xi1 = xi2.-1;
slope = ( log(y1(xi2)).-log(y1(xi1)) )./( x(xi2).-x(xi1) );
plot(x(xi2), slope );

# calculate total probability, should sum to one
total_prob1 = sum(y1)*1/xpart*xrng
total_prob2 = sum(y2)*1/xpart*xrng
sum(y2)*1/xpart*xrng - total_prob1
sum(y3)*1/xpart*xrng - total_prob1

# check convergence

slog_dist_simpson(  0, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 10, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 20, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 40, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 60, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 80, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 100, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 200, ita, beta, gr, phi, 4, 10)
slog_dist_simpson( 300, ita, beta, gr, phi, 4, 10)
