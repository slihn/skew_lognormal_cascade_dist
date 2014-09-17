
function q = slog_sqp_fn_hess ( pm )

    ita = pm(1,1); beta = pm(2,1); mu = pm(3,1); astd = pm(4,1); 

    q = astd;
    if ( ita >= 0.99 || ita <= 0 )
      q = -1;
      return;
    endif
    if ( beta >= 1.9 || beta <= 0 )
      q = -1;
    endif
endfunction
