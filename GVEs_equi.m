function statedot = GVEs_equi(t, state, mu, perturb, m0, mdot)
% GVES_EQUINOCTIAL  Calculates state derivative via equinoctial GVEs.
%
%   INPUT: 
%       t: current time step in integration, unused here [s]
%       state = [p f g h k L]': current state, [km, rad]
%       mu: gravitational parameter of central body [km^3/s^2]
%       perturb = [f_R f_T f_N]': perturbing force 
%           in RTN-frame, [kg*km/s^2]
%       m0: wet mass of satellite [kg]
%       mdot: rate of change of mass of sat during C.T. [kg/s]
%   OUTPUT:
%       statedot = d/dt (state), [km/s,rad/s]
    
    m = m0 - mdot * t;
    perturb = perturb / m;
    
    p = state(1); f = state(2); g = state(3);
    h = state(4); k = state(5); L = state(6);
    
    R = perturb(1); T = perturb(2); N = perturb(3);
    
    X = 1 + h^2 + k^2;
    w = 1 + f * cos(L) + g * sin(L);
    beta = sqrt(p/mu);
    
    statedot = zeros(6,1);
    
    statedot(1) = 2 * p * beta / w * T;
    
    statedot(2) = beta * (R*sin(L) + ((w+1)*cos(L)+f)/w*T...
        - g*(h*sin(L)-k*cos(L))/w*N);
    
    % The equations for gdot differ in the literature...proceed with caution.
    statedot(3) = beta * (-R*cos(L) + ((w+1)*sin(L)+g)/w*T...
        + f*(h*sin(L)-k*cos(L))/w*N); % JPL
    
    %statedot(3) = beta * (-R*cos(L) + ((w+1)*sin(L)+g)/w*T...
    %    + f*(h*sin(L)+k*cos(L))/w*N); % ESOC
    
    statedot(4) = beta * X/(2*w)*cos(L)*N;
    
    statedot(5) = beta * X/(2*w)*sin(L)*N;
    
    statedot(6) = sqrt(mu*p)*(w/p)^2 + beta*(h*sin(L)-k*cos(L))/w*N;
end