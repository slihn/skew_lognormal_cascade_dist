
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

# Solve zp and wd (sigma) for the skew case
# Some part of it uses brute force root-finding
# Not very elegant, but works

function [ zp, wd ] = slog_solve_zp_skew ( x, ita, beta, gr, phi ) 

debug = 0;

[ zp_1, wd_1 ] = slog_solve_zp_sym ( x, ita, phi );
if ( beta == 0.0 && gr == 0.0 ) 
    zp = zp_1;
    wd = wd_1;
    return;
endif

it_1 = slog_dist_integrand( zp_1, x, ita, beta, gr, phi );

# incorporate gr using iterative method

cnt=0;
while (1)
    hsym = zp_1-ita;
    fh_sym = phi * e.^(ita*hsym) .* gr;
    [ zp_2, wd_2 ] = slog_solve_zp_sym ( x-fh_sym, ita, phi );
    it_2 = slog_dist_integrand( zp_2, x, ita, beta, gr, phi );
    if ( it_2 > it_1 && abs(zp_1-zp_2)/wd_1 > 0.02 ) 
        zp_1 = zp_2;
        wd_1 = wd_2;
        it_1 = it_2;
        cnt++;
    else
        break;
    endif
endwhile

zp = zp_1;
wd = wd_1;
it = slog_dist_integrand( zp, x, ita, beta, gr, phi );
if ( debug )
    printf ("iterative cnt = %d: zp %f wd %f it %f\n", cnt, zp, wd, it);
endif

# incorporate beta

if ( gr <= 4 )
  prec_list = [ 0.1, 0.01, 0.001 ];
else
  prec_list = [ 0.4, 0.1, 0.01, 0.001 ];
endif

for prec = prec_list
    wd_change = 1;
    while (wd_change)

        delta = wd*prec;
        it_pos = slog_dist_integrand( zp+delta, x, ita, beta, gr, phi );

        if ( it_pos > it )
          direction = 1;
        else
          direction = -1;
        endif

        if ( debug )
                printf("precision %f zp %f wd %f delta %f dir %d\n", prec, zp, wd, delta, direction);
        endif

        cnt=0;
        while (1)
            zp_2 = zp + direction*delta;
            it_2 = slog_dist_integrand( zp_2, x, ita, beta, gr, phi );
            if ( it_2 > it )
                zp = zp_2;
                it = it_2;
                cnt++;
            else
                break;
            endif
        endwhile

        if ( debug )
            printf ("movement cnt = %d: zp %f wd %f it %f\n", cnt, zp, wd, it);
        endif

        wd_change = 0;
        cnt = 0;
        while (prec <= 0.10)
            it_wd1 = slog_dist_integrand( zp+wd, x, ita, beta, gr, phi ) / it;
            it_wd2 = slog_dist_integrand( zp-wd, x, ita, beta, gr, phi ) / it;
            itm = max([it_wd1,it_wd2]);
            if ( itm < 0.60 )
                if ( itm < 0.01 )
                    wd_change = 0.50;
                elseif ( itm < 0.40 )
                    wd_change = 0.80;
                else
                    wd_change = 0.90;
                endif
                wd *= wd_change;
                cnt++;
                if ( debug )
                    printf ("wd change cnt = %d: zp %f wd %f itmax %f\n", cnt, zp, wd, itm);
                endif
            else
                break;
            endif
        endwhile

    # if wd is shrinked, goto to delta=... again
    endwhile
endfor

if ( debug )
    printf("final: zp %f wd %f it %f\n", zp, wd, it);
endif

endfunction

