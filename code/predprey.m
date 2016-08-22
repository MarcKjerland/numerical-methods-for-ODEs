% Predator-Prey system
% also known as the Lotka-Volterra equations

function fu = predprey(u)

    % For use in 2d case only!
    if (length(u) ~= 2)
        'Dimension not equal to 2!'
        fu = u.*0; return
    end

    % u(1) is the prey population
    % u(2) is the predator population

    % Parameters
    a = 4.0; % growth rate of prey (food is always available)
    b = 2.0; % death rate of prey (as eaten by predators)
    c = 1.0; % growth rate of predators (their food is the prey)
    d = 1.0; % death rate of predators (from natural causes)

    fu = zeros(2,1);
    fu(1) = (a - b*u(2)) * u(1);
    fu(2) = (c*u(1) - d) * u(2);

end%function

