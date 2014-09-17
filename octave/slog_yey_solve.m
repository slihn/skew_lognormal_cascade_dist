
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

# solving z0=e^(-a0*z0)
# input: a0 can be an array, debug=1 if you want to see progress printed out
# output: [ z0, er, cnt ] or z0 only
#         z0 is the solution, er is the error, cnt is number of iterations

function varargout = slog_yey_solve(a0, debug)

  [ p,q ] = size(a0);

  z0 = zeros(p,q);
  er = zeros(p,q);
  cnt = zeros(p,q);
  for p2 = 1:p
  for q2 = 1:q
    [ z0(p2,q2), er(p2,q2), cnt(p2,q2) ] = slog_yey_solve1 ( a0(p2,q2), debug ); 
  endfor
  endfor

  if ( nargout == 1 )
    varargout{1} = z0;
    return;
  endif
  if ( nargout == 3 )
    varargout{1} = z0;
    varargout{2} = er;
    varargout{3} = cnt;
    return;
  endif
  error("undefined output pattern\n");

endfunction

  