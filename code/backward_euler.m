% This is the Backward Euler method
%  for a time-independent RHS.

% Given the RHS of a dynamical system f and the initial state u
%  i.e. u'(t) = f(u), 
%  this function computes one step of the B.E. method
%  and returns the update to u.

function u = backward_euler(f,u,dt)

% Want to solve y = u + dt*f(y) for y.
% Ideally we would use Newton's method
% to find the zeros of 
% ff(y) = -y + u + dt*f(y).
ff = @(y) -y + u + dt*f(y);

% Since the Jacobian df/dy isn't known,
% a quasi-Newton method must be used.

% Here is using Broyden's method:
tol=0.0001; % convergence tolerance
u = broyden(u,ff,tol);

% Matlab's built-in solver for
% nonlinear algebraic equations:
%u = fsolve(ff,u);

end %backward_euler
