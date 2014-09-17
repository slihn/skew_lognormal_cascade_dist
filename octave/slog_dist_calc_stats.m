
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

function f2_st = slog_dist_calc_stats ( ita, beta, gr, phi ) 
    
    f2_st = zeros(1,1);
    r  = beta*beta*ita*ita + 1.0;
    e1 = e^(ita*ita/2.0);
    e2 = e^(ita*ita*2.0);
    e3 = e^(ita*ita*9.0/2.0);
    e4 = e^(ita*ita*8.0);
    b1 = beta*ita*ita+gr;
    b2 = beta*ita*ita*2.0+gr;
    b3 = beta*ita*ita*3.0+gr;
    b4 = beta*ita*ita*4.0+gr;
    u1 = phi* e1*b1;
    k1 = u1;
    u2 = phi^2* e2*(r+b2*b2);
    u3 = phi^3* e3*b3* (3.0*r+b3*b3);
    u4 = phi^4* e4* (3.0*r*r+6.0*b4*b4*r+b4*b4*b4*b4);
    f2_st(1) = u1;
    k2 = u2-u1*u1;
    f2_st(2) = sqrt(k2); # std = sqrt(variance)
    k3 = u3-3.0*k2*k1-k1*k1*k1;
    f2_st(3) = k3/sqrt(k2*k2*k2);
    k4 = u4-4.0*k3*k1-3.0*k2*k2-6.0*k2*k1*k1-k1*k1*k1*k1;
    f2_st(4) = k4/k2/k2;

endfunction

