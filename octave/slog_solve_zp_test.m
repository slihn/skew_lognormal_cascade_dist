
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

more off;

phi=  1;
x    = 4
ita  = 0.9
beta = -0.9
gr   = 5
printf "----- begin ----------------------\n";

[ zp, wd ] = slog_solve_zp_skew ( x, ita, beta, gr, phi ) 

# brute force validation

n = 100000;
zn = (1:n);
zm1 = (1:(n-1));
zm2 = (2:n);

z  = (zn-n/2)/(n/10);
zd = (zm1-n/2)/(n/10);

h  = z-ita;
fh = phi * e.^(ita*h) .* (beta*ita*h+gr);

ft  = slog_dist_integrand( z, x, ita, beta, gr, phi );
ft2 = slog_dist_integrand( z, x, ita,    0, gr, phi );
ft3 = slog_dist_integrand( z, x, ita,    0,  0, phi );

fd = (ft(zm2)-ft(zm1))*1000; # 1:(n-1)

# visual check

figure(1);
subplot(2,1,1);
plot(z,ft/max(ft), z,ft2/max(ft2), z,ft3/max(ft3));
subplot(2,1,2);
plot(zd, fd);

# brute force root-finding

izero = -1;
i0 = min(find(fd==max(fd)));
for i = i0:(n-2)
  if ( fd(i+1) * fd(i) <= 0 )
    izero = i;
    break;
  endif
endfor

zp2 = zd(izero)

# error

err_zp = (zp2-zp)/zp
err_wd = (zp2-zp)/wd
