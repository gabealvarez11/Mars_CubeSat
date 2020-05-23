function statedot = GVEs_kep_esoc(t, state, mu, perturb, m0, mdot)
% GVES_KEP_ESOC  Calculates state derivative via ESOC's Keplerian GVEs.
%
%   INPUT: 
%       t: current time step in integration [s]
%       state = [a e i raan arg_periapsis nu]': current state, [km, rad]
%       mu: gravitational parameter of central body [km^3/s^2]
%       perturb = [f_R f_T f_N]': perturbing force 
%           in RTN-frame, [kg*km/s^2]
%       m0: wet mass of satellite [kg]
%       mdot: rate of change of mass of sat during C.T. [kg/s]
%   OUTPUT:
%       statedot = d/dt (state), [km/s,rad/s]
    
    m = m0 - mdot * t;
    perturb = perturb / m;
    
    a = state(1); e = state(2); i = state(3);
    raan = state(4); arg_periapsis = state(5); nu = state(6);
      
    R = perturb(1); T = perturb(2); N = perturb(3);
    
    n = sqrt(mu/a^3);
    p = a*(1-e^2);
    beta = sqrt(p/mu);
    r = p/(1+e*cos(nu));
    u = arg_periapsis + nu;
    h = sqrt(mu*p);
    
    statedot = zeros(6,1);
        
    statedot(1) = 2*a^2/h *(e*sin(nu)*R+(1+e*cos(nu))*T);

    statedot(2) = beta*(sin(nu)*R+(e+2*cos(nu)+e*cos(nu)^2)/(1+e*cos(nu))*T);
    
    statedot(3) = r/h*cos(u)*N;
    
    statedot(4) = r/h*sin(u)/sin(i)*N;
    
    statedot(5) = beta/e*(-cos(nu)*R+(2+e*cos(nu))/(1+e*cos(nu))*sin(nu)*T)...
        -cos(i)*statedot(4);
    
    statedot(6) = h/r^2 + beta/e*(cos(nu)*R-(2+e*cos(nu))/(1+e*cos(nu))*sin(nu)*T);
end