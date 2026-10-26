# ode67: Brainrot Numerical Integration

<p align="center">
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/651dfbd4-f06f-45ed-8657-3b961f8252e4" />
  <figcaption align="center">(I apologize for the AI generated image :sob:, I needed something incredibly brainrotted.)</figcaption>
</p>

`ode67` is a MATLAB implementation of Fehlberg's embedded Runge-Kutta 6(7) method from NASA TR R-287 [[1]](#references), with adaptive step-size control. Each step computes two solutions from the same ten stage evaluations: a 6th-order solution, which is used at the next time step, and a 7th-order solution that is used only to estimate the local truncation error. The step is accepted if the estimated local truncation error (difference between the two solutions) is within a user-set tolerance (or the default value), and the next step size is computed from that error estimate. The 6th-order propagated solution gives a global error that scales as $O(h^6)$.

The following Butcher Tableau was used (Page 48, [[1]](#references)):
<p align="center">
  <img width="600" alt="image" src="https://github.com/user-attachments/assets/77330df1-af0d-4e9d-b80f-e41b9992ba44" />
</p>

## Usage
After installing the script, make sure the folder it is in is added to your current workspace.

Usage of `ode67` is very similar to MATLAB's `ode45`. An example use case is shown here:
```matlab
rate_func = @(t,X) [X(2); -2*X(2) - 5*X(1)];
X0 = [1; 0];
tspan = [0 8];
max_local_err = 1e-6;
[t, X] = ode67(rate_func, tspan, X0, max_local_err);
```
The inputs are:
- `rate_func`: the system's rate function structured as a MATLAB anonymous function
- `tspan`: two-element vector representing the time interval of integration
- `X0`: the vector of initial conditions
- `max_local_err`: scalar defining the maximum allowed local truncation error (optional)

The outputs are:
- `t`: a vector of the resulting integration time steps
- `X`: numerical solution of all system states over the integration interval

Here is a plot result of `ode67` applied to the damped 2nd-order system from the code snippet above:
<p align="center">
  <img width="700" alt="demo_plot" src="https://github.com/user-attachments/assets/d59cbf67-f566-4b69-8608-a8410a93b9e5" />
</p>

## Analysis
A brief analysis on the global truncation error was performed:
<p align="center">
  <img width="700" alt="global_error" src="https://github.com/user-attachments/assets/91b34628-6943-4b27-ae31-14b82bb7a294" />
</p>

## References
1. E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control," NASA Technical Report TR R-287, 1968. [NTRS 19680027281](https://ntrs.nasa.gov/citations/19680027281).
