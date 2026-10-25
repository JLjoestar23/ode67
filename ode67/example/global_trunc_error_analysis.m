function global_trunc_error_analysis(rate_func, analytical_soln, tspan, num_trials)

    tol_list = logspace(-14, -8, num_trials);

    % storage for trial data
    h_avg_list    = zeros(1, num_trials);
    global_errors = zeros(1, num_trials);

    for i = 1:num_trials
        [global_errors(i), h_avg_list(i)] = calc_global_trunc_error(rate_func, analytical_soln, tspan, tol_list(i));
    end

    % central range for line fitting
    idx = max(1, floor(0.2*num_trials)) : min(num_trials, floor(0.8*num_trials));

    % error vs. average step size
    log_h   = log10(h_avg_list(idx));
    log_err = log10(global_errors(idx));
    pf      = polyfit(log_h, log_err, 1); % log(err) = slope*log(h) + b
    fit_h   = polyval(pf, log_h);

    figure;
    loglog(h_avg_list, global_errors, 'b.', 'MarkerSize', 15, ...
           'DisplayName', 'RK6(7) data');
    hold on;
    loglog(10.^log_h, 10.^fit_h, 'r--', 'LineWidth', 4, ...
           'DisplayName', sprintf('RK6(7) fit (slope = %.2f)', pf(1)));
    title('Global Error vs. Average Step Size');
    xlabel('Average timestep (s)');
    ylabel('Global error at end time');
    legend('Location', 'southeast');
    grid on; box on; hold off;
end