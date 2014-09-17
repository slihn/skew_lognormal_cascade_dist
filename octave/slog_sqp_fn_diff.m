
    # 1: sp: q1, and fit q2, skewness but not kurtosis
    # 2: dj: q1, fit skewness and kurtosis, but not q2
    # else: fit q1 only

function q = slog_sqp_fn_diff ( pm )
    global hx hy adlt hd ret_st fn_diff_mode hmin hmax;
    
    ita = pm(1,1); beta = pm(2,1); gr = pm(3,1); hfac = pm(4,1); delta = pm(5,1);
    
    printf ("at ita/beta = %.3f %.3f  gr/hfac %.3f %.3f dt %.3f\n", ita, beta, gr, hfac, delta);

    if ( ita >= 0.99 )
      q = abs(ita-0.99)*2000;
      return;
    endif
    if ( ita < 0.0 )
      q = abs(ita)*2000;
      return;
    endif
    #if ( beta >= 0.0 || beta <= -2.0 )
    #  q = abs(beta-1.0)*2000;
    #  return;
    #endif
    if ( abs(gr) > max(abs(hx)*hfac) )
      q = 1000;
      return;
    endif
    
    f2 = zeros(1,1);
    delta_max = hd(2)-hd(1);
    for xi = 1 : length(hd)
      f2(xi) = slog_dist_simpson((hd(xi)-delta)*hfac, ita, beta, gr, 1, 10, 10);
    endfor  
    f2 = f2 ./ sum(f2);
    f2 = log(f2);
    
    f1 = interp1(hx,hy,hd, 'spline', 'extrap');
    f1 = log(f1);
    zz = find(isnan(f1) .+ isinf(f1) .+ (f1<=-10));
    f1(zz) = -15;

    if ( length(find(isnan(f1))) > 0 )
      q = 1000;
      return;
      #error("f1 has NAN\n");
    endif    
    if ( length(find(isnan(f2))) > 0 )
      q = 1000;
      return;
      #error("f2 has NAN\n");
    endif    

    f2_st = slog_dist_calc_stats ( ita, beta, gr, 1/hfac ); 

    #if ( f2_st(1) < -2.0 )
    #  q = 10000*abs(u1);
    #  return;
    #endif
    
    figure(9);
    subplot(2,1,1);
    newplot;
    hold on;
    plot(hd*250, f1,".");
    plot(hx*250, log(hy),"r");
    plot(hd*250, f2,"b");
    title(sprintf("Q2,dbg dR %.2f ita/b %.3f %.3f %.3f %.3f mode %d",
          delta*250, ita,beta,gr,hfac, fn_diff_mode));
    hold off;

    subplot(2,1,2);
    newplot;
    hold on;
    plot(hd*250, e.^f1,".");
    plot(hx*250, hy,"y");
    plot(hd*250, e.^f2,"r");
    title(sprintf("Q1,dbg m %.3f %.3f std %.3f %.3f st8 %.2f %.2f st9 %.2f %.2f",
          ret_st(6)*250, f2_st(1)*250, 
          ret_st(7)*250, f2_st(2)*250, 
          ret_st(8),f2_st(3), ret_st(9),f2_st(4)));
    hold off;
    drawnow;

    printf("debug: m %.3f %.3f st8 %.2f %.2f st9 %.2f %.2f\n",
           ret_st(6)*250, f2_st(1)*250, ret_st(8),f2_st(3), ret_st(9),f2_st(4));
    
    f1(zz) = 0;
    f2(zz) = 0;
    
    zz2 = find( hd > hmax || hd < hmin );
    f1(zz2) = 0;
    f2(zz2) = 0;
    
    # normal region
    mid = find( abs(hd) <= ret_st(7) );
    
    # fit pdf
    q1 = 200*sqrt( sum((e.^f1-e.^f2).^2) / (length(hd)-length(zz)-1) );
    q1a = 1000*sqrt( sum((e.^f1(mid)-e.^f2(mid)).^2) / (length(hd)-1) );
    
    eh = e.^(hd-delta);
    q1e = 200*sqrt( sum((eh.*e.^f1-eh.*e.^f2).^2) / (length(hd)-length(zz)-1) );
    
    # fit log of pdf
    q2 = sqrt( sum((f1-f2).^2) / (length(hd)-length(zz)-1) );
    q2a = sqrt( sum((f1-f2).^2.*abs(hd)) / (length(hd)-length(zz)-1) );
    
    # skew and kurtosis
    # it was *40 and /5, but I changed for mkvalt fit 2/1/2009
    q3sk = abs( ret_st(8) - f2_st(3) ); 
    q3ku = abs( ret_st(9) - f2_st(4) ); 
    q3 = q3sk*3 + q3ku/5; 
    
    # mean and std
    # it is hard to give std a big weight, it messes up kurtosis
    q4 = abs(ret_st(6) - f2_st(1))*250 + abs(ret_st(7) - f2_st(2))*50;

    if ( fn_diff_mode == 1 )
      q = q1;
    elseif ( fn_diff_mode == 2 )
      q = q2;
    elseif ( fn_diff_mode == 3 )
      q = q1 + q2a;
    elseif ( fn_diff_mode == 4 ) # for mktcap fits
      q = q1 + q1e +q3 + q4;
    elseif ( fn_diff_mode == 5 ) # for low kurtosis
      q3 = q3sk*30 + q3ku*100;
      q = q1 + q3 + q4;
    elseif ( fn_diff_mode == 6 ) # for OU, std=2, kurtosis=6-10
      q = q1 + q2/10 + q4*10 + q3/10 + abs(ret_st(8));
    else
      q = q1 + q2 + q3 + q4;
    endif
    
    printf ("--%2d--------------- q = %16.10f %.4f %.4f %.5f %.5f\n", fn_diff_mode, q, q1,q2,q3,q4);

endfunction
