function [global_trunc_error, h_avg] = calc_global_trunc_error(rate_func_in, analytical_soln, tspan, max_local_err)
    
    % Get the current state from the analytical solution
    X0 = analytical_soln(tspan(1)); % start from exact initial condition
    
     % integrate with adaptive RK6(7); wrapper counts every rate_func call
    [t_list, X_list_num] = ode67(rate_func_in, tspan, X0, max_local_err);

    X_end_ana = analytical_soln(t_list(end));
    global_trunc_error = norm(X_list_num(:, end) - X_end_ana);
    h_avg = (tspan(2) - tspan(1)) / (numel(t_list) - 1);
end
