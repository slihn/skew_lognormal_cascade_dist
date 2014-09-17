
    # 1: sp: q1, and fit q2, skewness but not kurtosis
    # 2: dj: q1, fit skewness and kurtosis, but not q2
    # else: fit q1 only

function q = slog_sqp_so_fn_diff ( pm )

    global hx hy adlt hd ret_st fn_diff_mode hmin hmax;
    
    beta  = pm(1,1); gr    = pm(2,1); phi  = pm(3,1); 
    eta2  = pm(4,1); beta2 = pm(5,1); phi2 = pm(6,1); 
    delta = pm(7,1);
    
    annl = 250;
    
    printf ("at b/g/phi = %.3f %.3f %.3f 2)eta/b/phi %.3f %.3f %.3f dt %.3f\n", 
            beta, gr, phi*100, eta2, beta2, phi2, delta);

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
    
    f2 = zeros(1,1);
    delta_max = hd(2)-hd(1);
    for xi = 1 : length(hd(1,:))
      f2(xi) = slog_dist_so_simpson((hd(1,xi)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    endfor  
    f2 = f2 ./ sum(f2);
    f2 = log(f2);
    
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

    f2_st = slog_dist_so_calc_stats ( beta, gr, phi, eta2, beta2, phi2, 10, 20 ); 
    
    figure(9);
    subplot(2,1,1);
    newplot;
    hold on;
    plot(hd*annl, f1,".");
    plot(hx*annl, log(hy),"r");
    plot(hd*annl, f2,"b");
    title(sprintf("debug: dR %.2f b/g/phi = %.3f %.3f %.3f 2)eta/b/phi %.3f %.3f %.3f mode %d",
          delta*annl, beta, gr, phi*100, eta2, beta2, phi2, fn_diff_mode));
    hold off;

    subplot(2,1,2);
    newplot;
    hold on;
    plot(hd*annl, e.^f1,".");
    plot(hx*annl, hy,"y");
    plot(hd*annl, e.^f2,"r");
    axis([-f2_st(2)*10*annl, f2_st(2)*10*annl]);
    title(sprintf("debug: m %.3f %.3f std %.3f %.3f st8 %.2f %.2f st9 %.2f %.2f",
          ret_st(6)*annl, f2_st(1)*annl, 
          ret_st(7)*annl, f2_st(2)*annl, 
          ret_st(8),      f2_st(3), 
          ret_st(9),      f2_st(4)      ));
    hold off;
    drawnow;

    printf("debug: m %.3f %.3f s %.3f %.3f st8 %.2f %.2f st9 %.2f %.2f\n",
           ret_st(6)*annl, f2_st(1)*annl, 
           ret_st(7)*annl, f2_st(2)*annl, 
           ret_st(8),      f2_st(3), 
           ret_st(9),      f2_st(4)      );
    
    f1(zz) = 0;
    f2(zz) = 0;
    
    zz2 = find( hd > hmax || hd < hmin );
    f1(zz2) = 0;
    f2(zz2) = 0;
    
    # normal region
    mid = find( abs(hd) <= ret_st(7) );
    
    q1 = 200*sqrt( sum((e.^f1-e.^f2).^2) / (length(hd)-length(zz)-1) );
    q1a = 1000*sqrt( sum((e.^f1(mid)-e.^f2(mid)).^2) / (length(hd)-1) );
    
    q2 = 3* sqrt( sum((f1-f2).^2) / (length(hd)-length(zz)-1) );
    q2a = sqrt( sum((f1-f2).^2.*abs(hd)) / (length(hd)-length(zz)-1) );
    # it was *40 and /5, but I changed for mkvalt fit 2/1/2009
    q3 = abs( ret_st(8) - f2_st(3) )*40 + abs( ret_st(9) - f2_st(4) )/5; 
       
    # it is hard to give std a big weight, it messes up kurtosis
    q4 = abs(ret_st(6) - f2_st(1))*annl*30 + abs(ret_st(7) - f2_st(2))*annl*30;

    if ( fn_diff_mode == 1 )
      q = q1;
    elseif ( fn_diff_mode == 2 )
      q = q2;
    elseif ( fn_diff_mode == 3 )
      q = q1 + q2a;
    else
      q = q1 + q2 + q3 + q4;
    endif
    
    printf ("--%2d---- q = %16.10f %.4f %.4f %f %f\n", fn_diff_mode, q, q1,q2,q3,q4);

endfunction
