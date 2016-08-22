% Lorenz equations for use by ODE solver
%
% The Lorenz equations (1963) are given by:
%  x' = sigma*(y-x)
%  y' = x*(rho - z) - y
%  z' = x*y - beta*z
%
function dx = lorenz63(u)
    % For use in 3d only!
    if (length(u) ~= 3)
        'Dimension not equal to 3!'
        fu = u.*0; return
    end

    % These parameters will produce a chaotic system
    rho = 28; sigma = 10; beta = 8/3;
    dx = zeros(3,1);
    dx(1) = sigma*(u(2) - u(1));
    dx(2) = u(1)*(rho - u(3)) - u(2);
    dx(3) = u(1)*u(2) - beta*u(3);
    return
end

