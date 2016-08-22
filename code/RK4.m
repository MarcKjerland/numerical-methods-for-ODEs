% This is the standard Runge-Kutta 4th order scheme
%  for a time-independent RHS.
% Given the RHS of a dynamical system f and the initial state u
%  i.e. u'(t) = f(u), 
%  this function computes one step of RK4
%  and returns the update to u

function u = RK4(f,u,dt)

k1 = dt * f(u);
k2 = dt * f(u + k1/2);
k3 = dt * f(u + k2/2);
k4 = dt * f(u + k3);
    
u = u + (k1 + 2*k2 + 2*k3 + k4)/6;

end %function
