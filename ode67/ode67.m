% Adaptive timestep numerical integration using an embedded RK6(7) method.
% 
% INPUTS:
% rate_func_in: function used to compute dXdt
% tspan: two element vector that denotes the integration endpoints
% X0: initial condition vector
% p: how error scales with step size (error = k*h^p)
% error_desired: the desired local truncation error at each step
%
% OUTPUTS:
% t_list: time vector
% X_list: state vector
function [t_list, X_list] = ode67(rate_func, tspan, X0, max_local_err, p)

    arguments
        rate_func (1,1) function_handle
        tspan (1,:) double
        X0 (:,1) double
        max_local_err (1,1) double = 1e-6
        p (1,1) double = 7
    end
    
    % intialize
    t = tspan(1);
    t_end = tspan(2);
    t_list(1) = t;
    XA = X0(:);
    X_list = XA;
    h = (t_end - t) / 10; % initial timestep based on time interval
    
    % integrate within tspan
    while t < t_end
        
        % generate timestep to hit t_end exactly
        if t + h > t_end
            h = t_end - t;
        end

        % evaluate next step
        [XB, h_next, redo] = ode67_step(rate_func, t, XA, h, max_local_err, p);
        
        % if |XB1-XB2| > desired_error
        if redo
            % update timestep
            h = h_next;
            continue
        else
            t = t + h; % update time
            XA = XB; % update XA value
            X_list(:, end+1) = XB; % store result as row
            t_list(end+1) = t; % append new timestep
            
            % update timestep
            h = h_next;
        end
    end
end

% This function computes X(t+h) and the size of the succeeding time step. 
% It also determines whether to accept/reject the projected value of X(t+h)
% based on the estimated local truncation error.
function [XB, h_next, redo, XB1, XB2] = ode67_step(rate_func_in, t, XA, h, max_local_err, p)

    arguments
        rate_func_in (1,1) function_handle
        t (1,1) double
        XA (:,1) double
        h (1,1) double
        max_local_err (1,1) double = 1e-6
        p (1,1) double = 4
    end
    
    % Butcher Tableau for the 6(7) method
    C = [0, 2/33, 4/33, 2/11, 1/2, 2/3, 6/7, 1, 0, 1];
    B = [77/1440, 0, 0, 1771561/6289920, 32/105, 243/2560, 16807/74880, 11/270, 0, 0;
         11/864, 0, 0, 1771561/6289920, 32/105, 243/2560, 16807/74880, 0, 11/270, 11/270];
    A = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         2/33, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 4/33, 0, 0, 0, 0, 0, 0, 0, 0;
         1/22, 0, 3/22, 0, 0, 0, 0, 0, 0, 0;
         43/64, 0, -165/64, 77/32, 0, 0, 0, 0, 0, 0;
         -2383/486, 0, 1067/54, -26312/1701, 2176/1701, 0, 0, 0, 0, 0;
         10077/4802, 0, -5643/686, 116259/16807, -6240/16807, 1053/2401, 0, 0, 0, 0;
         -733/176, 0, 141/8, -335763/23296, 216/77, -4617/2816, 7203/9152, 0, 0, 0;
         15/352, 0, 0, -5445/46592, 18/77, -1215/5632, 1029/18304, 0, 0, 0;
         -1833/352, 0, 141/8, -51237/3584, 18/7, -729/512, 1029/1408, 0, 1, 0];

    % # of stages = order of the method
    stages = length(C);

    % state dimension
    n = length(XA);

    % pre-allocate vector of approximations
    % length of vector should match the order
    K = zeros(n, stages);
    
    % compute all K_i stages
    for i=1:stages
        % evaluate the sum of a_{i,j}*k_j terms as dot product
        sum_val = K*(A(i, :)');
        % evaluate the ith K approx
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*(sum_val));
    end
    
    % compute both embedded next timestep estimates
    XB1 = XA + h * (K*B(1,:)');
    XB2 = XA + h * (K*B(2,:)');

    % various constants
    alpha = 6; % step size increase factor cap
    
    % approximate the local truncation error
    err = norm(XB1-XB2);
    
    % update timestep
    h_next = min(alpha, 0.9*(max_local_err/err).^(1/p)) * h;
    
    if err > max_local_err
        XB = XA; % don't update next step estimate
        redo = true; % trigger a recalculation
    else
        % return XB if h is sufficient
        XB = XA + h * (K*B(1,:)');
        redo = false;
    end

end