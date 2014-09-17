
function pm_out = slog_fit_and_plot( ret_dump, h_bars, slog_param, fit_param, figure_no, axis_arr, text_pos, ttl )
    
    if ( length(slog_param) != 6 )
        error("slog_param length is not 6");
    endif
    if ( length(fit_param) != 4 )
        error("fit_param length is not 4");
    endif
    # ---------------------------------------------------------
    
    eta   = slog_param(1); beta = slog_param(2);
    gr    = slog_param(3); hfac = slog_param(4);
    delta = slog_param(5); annl = slog_param(6);

    do_sqp_fit     = fit_param(1);
    fn_diff_mode   = fit_param(2);
    R_cutoff_min   = fit_param(3);
    R_cutoff_max   = fit_param(4);
      
    if ( do_sqp_fit != 0 )  
        [ eta, beta, gr, hfac, delta ] = slog_sqp_standard ( ret_dump, h_bars, 
            eta, beta, gr, hfac, delta, 
            R_cutoff_min/annl, R_cutoff_max/annl, fn_diff_mode )
    endif
    # ---------------------------------------------------------

    st = statistics(ret_dump*annl)
      
    h = make_dist(ret_dump, h_bars);
    h_delta = ( max(h(1,:))-min(h(1,:)) ) / (h_bars-1);
    hi = find(h(2,:)>0);

    f2 = zeros(1,1);
    for xi = 1 : length(h(1,:))
      f2(xi) = slog_dist_simpson((h(1,xi)-delta), eta, beta, gr, 1/hfac, 10, 10);
    endfor  
    f2 = f2 ./ sum(f2);
    f2_st = slog_dist_calc_stats ( eta, beta, gr, 1/hfac ); 
    
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
        sprintf("Fit:\neta=%.2f\nbeta=%.2f\ngr=%.3f\nphi=%.3f\nstd=%.3f\nData statistics:\nmean=%.3f\nstd=%.3f\nskewness=%.2f\nkurtosis=%.2f",
        eta, beta, gr, 1/hfac*annl, f2_st(2)*annl,
        st(6), st(7), st(8), st(9) ));
    hold off;

    pm_out = [ eta, beta, gr, hfac, delta, annl ];
    
endfunction
