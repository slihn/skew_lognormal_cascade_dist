
# hfac is 1/phi 

function [ ita, beta, gr, hfac, delta ] = slog_sqp_standard ( 
    arr, bars, 
    ita_in, beta_in, gr_in, hfac_in, delta_in, 
    hmin_in, hmax_in, md )

  global ret_st fn_diff_mode hmin hmax;
  
  fn_diff_mode = md;
  hmin = hmin_in;
  hmax = hmax_in;
  
  h = make_dist(arr, bars);
  # h = slog_sqp_trim_tail ( h, bars);
  [ mu, astd, h ] = hist_fit_normal ( h );
    
  cont_set_global_h ( h );
  ret_st = statistics( arr )

  i0 = [ ita_in; beta_in; gr_in; hfac_in; delta_in ]
  [fnl, obj, info, iter, nf, lambda] = sqp(i0, 
   #  @slog_sqp_fn_diff, @slog_sqp_fn_grad, @slog_sqp_fn_hess)
      @slog_sqp_fn_diff )
  
  ita   = fnl(1,1);
  beta  = fnl(2,1);
  gr    = fnl(3,1);
  hfac  = fnl(4,1);
  delta = fnl(5,1);
    
endfunction

function h = slog_sqp_trim_tail ( h, bars)

  n = length(h(1,:));
  mi = bars/2;
  mh = h(2,mi);
  for i = 1:2:(n-1)
    r = i:(i+1);
    f = (h(2,r) > 0);
    if ( sum(f) != 2 )
      h(2,i) = 0;
      h(2,i+1) = 0;
    endif
    if ( h(2,i) > mh )
      mi = i;
    endif
  endfor
  
  # trim twice
  if ( 1 )
    for i = 1:2:(n-1)
      r = i:(i+1);
      f = (h(2,r) > 0);
      if ( sum(f) != 2 )
        h(2,i) = 0;
        h(2,i+1) = 0;
      endif
    endfor
  endif
  
  ci = 1;
  for i = 1:mi
    if ( h(2,i) == 0 )
      ci = i;
    endif
  endfor
  for i = mi:bars
    if ( (i-mi)/(bars-mi) > (mi-ci)/mi )
      h(2,i) = 0;
    endif
  endfor
  printf("tails trimmed...\n");
  pause(5);
endfunction
