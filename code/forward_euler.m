% This is the Forward Euler method
%  for a time-independent RHS.

% Given the RHS of a dynamical system f and the initial state u
%  i.e. u'(t) = f(u), 
%  this function computes one step of the F.E. method
%  and returns the update to u.

function u = forward_euler(f,u,dt)

u = u + dt * f(u);
    
end %function



