
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
# scalar input only

# ------------------------------------------------------  
function [ z0, er, cnt ] = slog_yey_solve1 (a0, debug)

  if ( a0 <= 0 )
    if ( a0 == 0 )
      z0 = 1;
      er = 0;
      cnt = 1;
      return;
    endif
    error("slog_yey_solve: a0 is negative");    
  endif

  amid = log(4);
  ag = log(a0/amid);
  d=0.83+1./(1.+e.^(ag/(0.8)))*0.15;
  fn = 1./(1.+e.^(d.*ag)); # approximate z
  fn2 = -0.07*ag.*e.^(-0.02*abs(ag).^2).*(0.5-abs(fn-0.5));
  fn  = fn.-fn2;
  
  
  if ( a0 < 0.1 )
    z0 = fn;
    cnt = 0;
    do
      z1 = e^(-a0*z0);
      er = abs((z1-z0)/z0);  
      cnt++;
      if ( debug )
        printf("prog1 %d %.4f\n",cnt,er);
      endif
      z0 = z1;
    until (er < 0.0001);

  elseif ( a0 > e^20 )  
    p = log(a0)-1;
    q0 = p-log(1+p-log(1+p));    
    y0 = e.^(1+q0);
    cnt = 0;
    do
      y1 = a0/log(y0);
      er = abs((y1-y0)/y0);  
      cnt++;
      if ( debug )
        printf("prog2 %d %.4f\n",cnt,er);
      endif
      y0 = y1;
    until (er < 0.0001);
    z0 = 1/y0;

  else
    cnt = 0;
    p = log(a0)-1;
    if ( p <= 1 )
      if ( debug )
        printf("prog3 quadratic init\n");
      endif
      b = (e^p+1)/e^p;
      c = 2*(e^p-1)/e^p;
      q0 = b-sqrt(b^2-c);
    elseif ( p <= 10 )
      if ( debug )
        printf("prog3 fn init\n");
      endif
      q0 = -log(fn)-1;
    else
      q0 = p-log(1+p-log(1+p));    
    endif

    z0 = e.^(-1.-q0);
    if ( abs(q0) > 0.0001 )
      do
        q1 = (e^(p-q0)*(1+q0)-1)/(e^(p-q0)+1);
        er = abs((q1-q0)/q0);  
        cnt++;
        if ( debug )
          printf("prog3 %d %.4f\n",cnt,er);
        endif
        q0 = q1;
      until (er < 0.0001);
    endif
    z0 = e.^(-1.-q0);

  endif
  
  # verify
  z1 = e^(-a0*z0);
  er = abs((z1-z0)/z0);  
  if ( er > 0.0001 )
    error(sprintf("slog_yey_solve: error too high: z=%f on a=%f err %f\n", z0,a0,er));
  endif
  if ( debug )
    printf("-end- %d %.4f (final error estimate in z)\n",cnt,er);
    a0
    z0
  endif
endfunction

  