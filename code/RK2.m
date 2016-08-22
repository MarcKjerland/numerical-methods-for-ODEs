% This is the midpoint method
%  also known as RK2 (it's a 2nd-order Runge-Kutta scheme)
%  for a time-independent RHS.
% Given the RHS of a dynamical system f and the initial state u
%  i.e. u'(t) = f(u), 
%  this function computes one step of the midpoint method
%  and returns the update to u.

function u = RK2(f,u,dt)

k1 = dt * f(u);
k2 = dt * f(u + k1/2);
    
u = u + k2;

% (the unusual choice of variables is meant for 
%   direct comparison with RK4)

end %function

