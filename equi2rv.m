function inertial_rv = equi2rv(state, mu)
% EQUI2RV  Converts equinoctial elements to intertial position, velocity.
%   
%   INPUT: equi = [p f g h k L], in [km, rad]
%          mu [km^3/s^2]
%   OUTPUT: inertial_rv = [r v]', in [km, km/s]

    p = state(1); f = state(2); g = state(3);
    h = state(4); k = state(5); L = state(6);
    
    inertial_rv = zeros(6,1);
    
    alpha2 = h^2 - k^2;
    s2 = 1 + h^2 + k^2;
    w = 1 + f * cos(L) + g * sin(L);
    r = p / w;
    gamma = sqrt(mu/p);
    
    inertial_rv(1) = r/s2 * (cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
    inertial_rv(2) = r/s2 * (sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
    inertial_rv(3) = 2*r/s2 * (h*sin(L) - k*cos(L));
    
    inertial_rv(4) = -gamma/s2 * (sin(L) + alpha2*sin(L)-2*h*k*cos(L)...
        + g - 2*f*h*k + alpha2*g);
    
    inertial_rv(5) = -gamma/s2 * (-cos(L) + alpha2*cos(L) + 2*h*k*sin(L)...
        - f + 2*g*h*k + alpha2*f);
    
    inertial_rv(6) = 2*gamma/s2 * (h*cos(L) + k*sin(L) + f*h + g*k);
end