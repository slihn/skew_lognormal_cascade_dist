
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

# n=10, nsigma=20 is good for DOW parameters, eta=0.75
#       slog_dist_so_calc_stats ( -0.1, 0.08, 1.6, 0.1, 0.0, 0.7, 10, 20)
# n=20, nsigma=40 is good for eta up to 1.0

function f2_st = slog_dist_so_calc_stats ( b, g, phi, eta2, beta2, phi2, n, nsigma ) 
    
    f2_st = zeros(1,1);

    r2  = beta2*beta2*eta2*eta2 + 1.0;
    eta = phi2*e^(eta2^2)*sqrt(r2);
    gr2 = -beta2 * eta2^2;

    # handle phi=0 i.e. eta=0
    if ( eta == 0 )
        f2_st = slog_dist_calc_stats ( eta, b, g, phi, n, nsigma );
        return;
    endif
    
    # eta is the std of the integrand, H = +/- eta* (20 or 40 if eta > 0.5)
    dH = eta / n;
    m  = n * nsigma;
    Hi = ( (-m) : 1 :  m );
    H  = Hi * dH;
  
    Hi_odd  = find(abs(rem(Hi,2)) == 1);
    Hi_even = find(rem(Hi,2) == 0);  

    p_H = slog_dist_simpson(H, eta2, beta2, gr2, phi2, 20, 20);

    # u0i = phi*p_H;
    # sum(u0i)*dH  # debug, this should be one...
    
    u1i = phi*(b.*H+g) .*e.^H .*p_H;
    u1  = dH/3 .* ( sum(u1i(Hi_odd))*2 + sum(u1i(Hi_even))*4 );

    u2i = phi^2*(b^2*H.^2.+2*b*g.*H.+g^2.+1) .*e.^(2*H) .*p_H;
    u2  = dH/3 .* ( sum(u2i(Hi_odd))*2 + sum(u2i(Hi_even))*4 );
    
    u3i = phi^3*(b.*H+g) .*(b^2.*H.^2.+2*b*g.*H.+g^2.+3) .*e.^(3*H) .*p_H;
    u3  = dH/3 .* ( sum(u3i(Hi_odd))*2 + sum(u3i(Hi_even))*4 );;
    
    u4i = phi^4*(b^4.*H.^4.+4*b^3*g.*H.^3.+6*b^2*g^2.*H.^2.+6*b^2.*H.^2.+4*b*g^3.*H.+12*b*g.*H.+g^4.+6*g^2.+3) .*e.^(4*H) .*p_H;
    u4  = dH/3 .* ( sum(u4i(Hi_odd))*2 + sum(u4i(Hi_even))*4 );
    
    k1 = u1;
    f2_st(1) = u1;
    k2 = u2-u1*u1;
    f2_st(2) = sqrt(k2); # std = sqrt(variance)
    k3 = u3-3.0*k2*k1-k1*k1*k1;
    f2_st(3) = k3/sqrt(k2*k2*k2);
    k4 = u4-4.0*k3*k1-3.0*k2*k2-6.0*k2*k1*k1-k1*k1*k1*k1;
    f2_st(4) = k4/k2/k2;
    
endfunction

