
    # 1: sp: q1, and fit q2, skewness but not kurtosis
    # 2: dj: q1, fit skewness and kurtosis, but not q2
    # else: fit q1 only

function q = slog_sqp_so_fn_diff2 ( pm )

    global hx hy adlt hd ret_st fn_diff_mode hmin hmax;
    
    beta  = pm(1,1); gr    = pm(2,1); phi  = pm(3,1); 
    eta2  = pm(4,1); beta2 = pm(5,1); phi2 = pm(6,1); 
    delta = pm(7,1);
    
    annl = 250;
    
    printf ("at b/g/phX = %.3f %.3f %.3f/%.2f [2]e/b/ph %.3f %.3f %.3f dR %.3f\n", 
            beta, gr, phi*100, 1/phi, eta2, beta2, phi2, delta*annl);

    eta = phi2*e^(eta2^2);

    if ( eta >= 0.99 || eta < 0.0 )
      q = abs(eta-0.5)*2000;
      return;
    endif
    if ( eta2 >= 0.5 || eta2 < 0.0 )
      q = abs(eta2-0.5)*2000;
      return;
    endif
    if ( beta2 >= 0.8 || beta2 < -0.8 )
      q = abs(beta2)*2000;
      return;
    endif

    if ( abs(gr) > max(abs(hx)/phi) )
      q = 1000;
      return;
    endif
    
    xi_max = find(hy==max(hy));
    max_h = hy(xi_max);

    h_delta = hd(1,2)-hd(1,1);
    max_f2 = h_delta*slog_dist_so_simpson((hd(1,xi_max)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    printf ("max_pair = %d %.4f f2 %.4f\n", xi_max, max_h, max_f2);

    hy1p = hy(xi_max+1); 
    hy1n = hy(xi_max-1);
    hy2p = hy(xi_max+2);
    hy2n = hy(xi_max-2);
    f2_1p = h_delta*slog_dist_so_simpson((hd(1,xi_max+1)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    f2_1n = h_delta*slog_dist_so_simpson((hd(1,xi_max-1)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    f2_2p = h_delta*slog_dist_so_simpson((hd(1,xi_max+2)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    f2_2n = h_delta*slog_dist_so_simpson((hd(1,xi_max-2)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    
    f2_st = slog_dist_so_calc_stats ( beta, gr, phi, eta2, beta2, phi2, 10, 20 ); 
    
    printf("debug: dR %.2f b/g/phX = %.3f %.3f %.3f 2)eta/b/phi %.3f %.3f %.3f mode %d\n",
          delta*annl, beta, gr, phi*100, eta2, beta2, phi2, fn_diff_mode);

    printf("debug: m %.3f %.3f std %.3f %.3f st8 %.2f %.2f st9 %.2f %.2f\n",
           ret_st(6)*annl, f2_st(1)*annl, 
           ret_st(7)*annl, f2_st(2)*annl, 
           ret_st(8),      f2_st(3), 
           ret_st(9),      f2_st(4)      );
        
    # normal region
    mid = find( abs(hd) <= ret_st(7) );
    
    q3 = abs( ret_st(8) - f2_st(3) )*40 + abs( ret_st(9) - f2_st(4) ) +abs(beta2)/1000; 
       
    # it is hard to give std a big weight, it messes up kurtosis
    q4 = abs(ret_st(6) - f2_st(1))*annl*100 + abs(ret_st(7) - f2_st(2))*annl*30;
    q5 = 50*(abs(max_h-max_f2 ) +abs(hy1p-f2_1p) +abs(hy1n-f2_1n) +abs(hy2p-f2_2p) +abs(hy2n-f2_2n) );
    
    q = q3 + q4 + q5;
    
    printf ("--%2d---- q = %16.10f %.4f %.4f %.4f\n\n", fn_diff_mode, q, q3,q4,q5);

endfunction
