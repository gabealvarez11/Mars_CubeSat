% Accepts as inputs eccentricity e and true anomaly at epoch nu in radians.
% Returns mean anomaly M in radians.
function M = nu2M(e, nu)
    % Find the eccentric anomaly E.
    cosE = (e + cos(nu)) / (1 + e * cos(nu));
    sinE = sin(nu) * sqrt(1 - e^2) / (1 + e * cos(nu));

    if sinE > 0
        E = acos(cosE);
    else % sinE < 0
        E = 2 * pi - acos(cosE);
    end
    
    % Apply the Kepler Equation.
    M = E - e * sin(E);
end