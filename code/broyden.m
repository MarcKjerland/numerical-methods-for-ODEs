% Broyden's method for solving
%  f(x) = 0,
% a system of n nonlinear equations in n variables.

% Requires an initial guess for the solution x,
% a user-defined function f,
% and accuracy tolerance tol.
% Returns a solution vector xv.
% WARNING. Method may fail to converge.

% This code is adapted from John Penny, Numerical Methods using Matlab, 2e.
% https://www.mathworks.com/matlabcentral/fileexchange/2305-numerical-methods-using-matlab-2e/content/edition2/na_funcs/broyden.m

function xv=broyden(x,f,tol)
xv=x;
fr=f(xv);
%Set initial inverse Jacobian
n=length(xv);
Jinv=eye(n);
it=0;
while norm(fr)>tol
    it=it+1;
    % One iteration of Newton's method
    pr=-Jinv*fr;
    xv=xv+pr;
    oldfr=fr;
    fr=f(xv);
    % Update approximation to inverse Jacobian
    % using Broydens formula
    y=fr-oldfr;
    oldJinv=Jinv;
    oyp=oldJinv*y-pr;
    pB=pr'*oldJinv;
    for i=1:n
        for j=1:n
            M(i,j)=oyp(i)*pB(j);
        end;
    end;
    Jinv=oldJinv-M./(pr'*oldJinv*y);
    %Exit condition
    %in case it doesn't converge
    if (it > 10000)
        break
    end
end;

end%broyden

