
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

# integration using Simpson's Rule near the peak position (zp)

# n: intervals within a sigma, typically 5-8
# nsigma: how many simga (wd) to integrate. 8 is pretty wide

# The range of integration is z = [ zp-nsigma*wd, zp+nsigma*wd ]

function v = slog_dist_so_simpson1 ( x, beta, gr, phi, eta2, beta2, phi2, n, nsigma ) 

    # calculate zp, wd (sigma)
    eta = phi2*e^(eta2^2);
    gr2 = -beta2 * eta2^2;

    # handle phi=0 i.e. eta=0
    if ( eta == 0 )
        v = slog_dist_simpson( x, eta, beta, gr, phi, n, nsigma );
        return;
    endif
  
    [ zp, wd ] = slog_solve_zp_skew ( x, eta, beta, gr, phi ); 
    h = wd / n;
    m = n * nsigma;
    zodd  = (   (-m) : 2 :  m    )*h + zp;
    zeven = ( (-m+1) : 2 : (m-1) )*h + zp;
  
    Hodd  = ( zodd  .- eta ) .* eta;
    Heven = ( zeven .- eta ) .* eta;

    ig2 = 2 * slog_dist_kernal( x,  Hodd, beta, gr, phi ) .* slog_dist_simpson(  Hodd, eta2, beta2, gr2, phi2, n, nsigma);
    ig3 = 4 * slog_dist_kernal( x, Heven, beta, gr, phi ) .* slog_dist_simpson( Heven, eta2, beta2, gr2, phi2, n, nsigma);
    v  = h/3 .* ( sum(ig2)+sum(ig3) ) * eta;

endfunction

