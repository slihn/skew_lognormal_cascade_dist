
function pm_out = slog_fit_and_plot_so( ret_dump, h_bars, slog_param, fit_param, figure_no, axis_arr, text_pos, ttl )
    
    if ( length(slog_param) != 8 )
        error("slog_param length is not 8");
    endif
    if ( length(fit_param) != 4 )
        error("fit_param length is not 4");
    endif
    # ---------------------------------------------------------
    
    beta  = slog_param(1); gr    = slog_param(2); phi  = slog_param(3);
    eta2  = slog_param(4); beta2 = slog_param(5); phi2 = slog_param(6);
    delta = slog_param(7); annl  = slog_param(8);

    do_sqp_fit     = fit_param(1);
    fn_diff_mode   = fit_param(2);
    R_cutoff_min   = fit_param(3);
    R_cutoff_max   = fit_param(4);
      
    if ( do_sqp_fit != 0 )  
        [ beta, gr, phi, eta2, beta2, phi2, delta ] = slog_sqp_so_standard ( ret_dump, h_bars, 
            beta, gr, phi, eta2, beta2, phi2, delta, 
            R_cutoff_min/annl, R_cutoff_max/annl, fn_diff_mode )
    endif
    # ---------------------------------------------------------
      
    h = make_dist(ret_dump, h_bars);
    h_delta = ( max(h(1,:))-min(h(1,:)) ) / (h_bars-1);
    hi = find(h(2,:)>0);

    f2 = zeros(1,1);
    for xi = 1 : length(h(1,:))
      f2(xi) = slog_dist_so_simpson((h(1,xi)-delta), beta, gr, phi, eta2, beta2, phi2, 10, 20);
    endfor  
    f2 = f2*h_delta;
    f2_prob = sum(f2)
    
    max_prob = [ max(h(2,:)), find(h(2,:)==max(h(2,:))), max(f2), find(f2==max(f2)) ]
    f2_st = slog_dist_so_calc_stats ( beta, gr, phi*annl, eta2, beta2, phi2, 10, 20 )
    st = statistics(ret_dump*annl);
    th_st = st(6:9)
    
    # ---------------------------------------------------------

    figure(figure_no);
    newplot();
    hold on;
    plot( h(1,:)*annl, log(h(2,:))-log(h_delta*annl), ".r");
    plot( h(1,:)*annl, log(f2)-log(h_delta*annl), "b");
    legend("data", "fit");
    xlabel(sprintf("R (x%d)",annl));
    ylabel("log(pdf)");
    title(ttl);
    if ( sum(abs(axis_arr)) != 0 )
        axis(axis_arr);  
    endif
    text(text_pos(1), text_pos(2), 
        sprintf("Fit:\nbeta=%.2f\ngr=%.3f\nphi=%.3f\neta2=%.3f\nbeta2=%.3f\nphi2=%.3f\nStatistics:\nmean=%.3f\nstd=%.3f\nskewness=%.2f\nkurtosis=%.2f",
        beta, gr, phi*annl, eta2, beta2, phi2,
        st(6), st(7), st(8), st(9) ));
    hold off;

    pm_out = [ beta, gr, phi, eta2, beta2, phi2, delta, annl ];
    
endfunction
