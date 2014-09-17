

function [ ita, beta, gr, hfac ] = slog_sqp_distfit ( ita_in, beta_in, gr_in, hfac_in, hmin_in, hmax_in, md )

  global fn_diff_mode hmin hmax;
  
  fn_diff_mode = md;
  hmin = hmin_in;
  hmax = hmax_in;
  
  i0 = [ ita_in; beta_in; gr_in; hfac_in ]
  [fnl, obj, info, iter, nf, lambda] = sqp(i0, 
      @slog_sqp_distfit_fn_diff )
  
  ita  = fnl(1,1);
  beta = fnl(2,1);
  gr   = fnl(3,1);
  hfac = fnl(4,1);
    
endfunction

