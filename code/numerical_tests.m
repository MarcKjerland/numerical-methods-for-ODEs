% This MATLAB script is a showcase of
% numerical integrators for ODEs
% for solving the initial value problem
%  u'(t) = f(u)
%  u(t0)  = u0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem Setup (edit this)
%
% Define the function on the right hand side (RHS)
% of the differential equation
%  (don't forget to change initial values)
%
% This function can be:
% 1. a Matlab built-in function
f = @sin;
% 2. a function defined in a function M-file
f = @predprey; % 2D
f = @lorenz63; % 3D
% 3. an anonymous function
f = @(z) -10*z; % exponential decay
f = @(z) [0 -1;1 0]*z; % sinusoid
% (In all cases, the function handle @ is needed)
%
% Let's use this function:
f = @(z) [0 -1;1 0]*z;

% Initial state (can be scalar or vector)
%u0 = 4; %scalar
u0 = [4;2]; % 2d vector
%u0 = [4;2;3]; %3d vector

% Timestep size
dt = 0.2; % important for stability!

% Start and stop times
t_start = 0;
t_stop = 10;

%
%%% end of problem setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we store the solution array u[i] = u(i*dt+t0)
% for each numerical integrator
% initialized as the starting state u0
ts = [t_start : dt : t_stop];
Nt = length(ts);
u_FE = u0*ones(1,Nt);
u_BE = u0*ones(1,Nt);
u_RK2 = u0*ones(1,Nt);
u_RK4 = u0*ones(1,Nt);

% Numerical integration
for i=1:Nt-1

    % Call explicit integrators
    u_FE(:,i+1) = forward_euler(f,u_FE(:,i),dt);
    u_RK2(:,i+1) = RK2(f,u_RK2(:,i),dt);
    u_RK4(:,i+1) = RK4(f,u_RK4(:,i),dt);

    % Call implicit integrators
    u_BE(:,i+1) = backward_euler(f,u_BE(:,i),dt);

end
figure
hold on
plot(ts,u_FE(1,:),'-')
plot(ts,u_BE(1,:),'--')
plot(ts,u_RK2(1,:),'-.')
plot(ts,u_RK4(1,:),':')
xlabel('t')
ylabel('U(t)')
legend('Forward Euler','Backward Euler','Midpoint method','Runge-Kutta 4','Location','Best')
title('Numerical methods comparison')
hold off

% Bonus: if you want to try the built-in RK45 solver
if (false)
    % ode45 expects a time-dependent system
    f45 = @(t,u) f(u);
    % The timestep is adaptive, so we only need
    % to specify the start and end times
    [t45 u45] = ode45(f45,[t_start t_stop],u0);
    figure; plot(t45,u45,'o')
    title('Time adaptive Runge-Kutta')
end
%

% Bonus 2: generates 3D plots i.e. for Lorenz 63 system
% For best results use t_stop=100 and dt=0.01
if (false)
    figure
    plot3(u_RK4(1,:),u_RK4(2,:),u_RK4(3,:))
    title('Lorenz 63 system')
end
