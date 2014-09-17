
    # 1: sp: q1, and fit q2, skewness but not kurtosis
    # 2: dj: q1, fit skewness and kurtosis, but not q2
    # else: fit q1 only

function q = slog_sqp_distfit_fn_diff ( pm )
    global hx hy adlt hd fn_diff_mode hmin hmax f2_st;
    
    ita = pm(1,1); beta = pm(2,1); gr = pm(3,1); hfac = pm(4,1); 
    
    printf ("at ita/beta = %f %f  gr/hfac %f %f\n", ita, beta, gr, hfac);

    if ( ita >= 0.99 || ita < 0.0 )
      q = abs(ita-0.5)*2000;
      return;
    endif
    #if ( beta >= 0.0 || beta <= -2.0 )
    #  q = abs(beta-1.0)*2000;
    #  return;
    #endif
    if ( abs(gr) > max(abs(hx)) )
      q = 1000;
      return;
    endif
    
    f2 = zeros(1,1);
    for xi = 1 : length(hd)
      f2(xi) = slog_dist_simpson(hd(xi)*hfac, ita, beta, gr, 1, 10, 10);
    endfor  
    f2 = f2 ./ sum(f2);
    f2 = log(f2)-log(0.1);
    
    f1 = interp1(hx,hy,hd, 'spline', 'extrap');
    f1 = log(f1);
    zz = find(isnan(f1) .+ isinf(f1) .+ (f1<=-10));
    f1(zz) = -15;

    if ( length(find(isnan(f1))) > 0 )
      error("f1 has NAN\n");
    endif    
    if ( length(find(isnan(f2))) > 0 )
      error("f2 has NAN\n");
    endif    

    f2_st = slog_dist_calc_stats ( ita, beta, gr, 1 );

    if ( f2_st(1) < -2.0 )
      q = 10000*abs(f2_st(1));
      return;
    endif
    
    figure(9);
    subplot(2,1,1);
    newplot;
    hold on;
    plot(hd, f1,".");
    plot(hx, log(hy),"r");
    plot(hd, f2,"b");
    title(sprintf("debug: ita/b %f %f %f %f mode %d",ita,beta,gr,hfac, fn_diff_mode));
    hold off;

    subplot(2,1,2);
    newplot;
    hold on;
    plot(hd, e.^f1,".");
    plot(hx, hy,"y");
    plot(hd, e.^f2,"r");
    title(sprintf("debug: st8 %.2f st9 %.2f",f2_st(3), f2_st(4)));
    hold off;
    drawnow;

    printf("debug: st8 %.2f st9 %.2f\n", f2_st(3), f2_st(4));
    
    f1(zz) = 0;
    f2(zz) = 0;
    
    zz2 = find( hd > hmax || hd < hmin );
    f1(zz2) = 0;
    f2(zz2) = 0;
    
    q1 = 70*sqrt( sum((e.^f1-e.^f2).^2) / (length(hd)-length(zz)-1) );
    q2 = sqrt( sum((f1-f2).^2) / (length(hd)-length(zz)-1) );
    # q3 = abs( ret_st(8) - f2_st(3) )*40 + abs( ret_st(9) - f2_st(4) )/5; 
    q3 = 0;
    
    if ( fn_diff_mode == 1 )
      q = q1;
    elseif ( fn_diff_mode == 2 )
      q = q2;
    else
      q = q1 + q2 + q3;
    endif
    
    printf ("-------------------------- q = %16.10f %f %f %f\n", q, q1,q2,q3);

endfunction
