% test ODE and its analytical solution
f1  = @(t,X) [X(2); -2*X(2) - 5*X(1)];
ex1 = @(t) [exp(-t).*(cos(2*t) + 0.5*sin(2*t)); ...
            -2.5*exp(-t).*sin(2*t)];

% setup ode67 test
X0    = ex1(0);
tspan = [0 8];
% max_local_err = 1e-6;
% p = 4;

[t, X] = ode67(f1, tspan, X0);

t_ana = linspace(tspan(1), tspan(2), 500);
X_ana = ex1(t_ana);
err = vecnorm(X - ex1(t), 2); % local truncation error

figure('Color', 'w');
subplot(2, 1, 1);
hold on;
plot(t_ana, X_ana(1,:), '-', 'Color', '#0b81a2', 'LineWidth', 2, 'DisplayName', 'x (analytical)');
plot(t_ana, X_ana(2,:), '-', 'Color', '#e25759', 'LineWidth', 2, 'DisplayName', 'x'' (analytical)');
plot(t, X(1,:), '.', 'Color', '#9d2c00', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'x (ode67)');
plot(t, X(2,:), '.', 'Color', '#003a7d', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'x'' (ode67)');
xlabel('Time (s)'); ylabel('State');
title(sprintf('Damped oscillator: RK6(7), %d timesteps', numel(t)-1));
legend('Location', 'northeast', 'NumColumns', 2);
grid on;
hold off;

subplot(2, 1, 2)
semilogy(t, max(err, eps), '.-', 'Color', '#0b81a2', 'LineWidth', 1.5, 'MarkerSize', 15);
xlabel('Time (s)'); ylabel('Local Truncation Error');
title('Local Truncation Error Per Timestep');
grid on; box off;


%% Global truncation error analysis

global_trunc_error_analysis(f1, ex1, [0 10], 500);