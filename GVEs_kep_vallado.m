function statedot = GVEs_kep_vallado(t, state, mu, perturb, m0, mdot)
% GVES_KEP  Calculates state derivative via Vallado's Keplerian GVEs (p. 636).
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
%
%   WARNING: Behaves strangely for along-track (T) perturbations.

    m = m0 - mdot * t;
    perturb = perturb / m;
    
    a = state(1); e = state(2); i = state(3);
    raan = state(4); arg_periapsis = state(5); nu = state(6);
      
    R = perturb(1); T = perturb(2); N = perturb(3);
    
    n = sqrt(mu/a^3);
    factor = sqrt(1-e^2);
    p = a*(1-e^2);
    r = p/(1+e*cos(nu));
    u = arg_periapsis + nu;
    h = sqrt(mu*p);
    
    statedot = zeros(6,1);
        
    statedot(1) = 2/(n*factor)*(e*sin(nu)*R + p/r*T);

    statedot(2) = factor/(n*a)*(sin(nu)*R + (cos(nu)+(e+cos(nu))/(1+e*cos(nu)))*T);

    statedot(3) = r*cos(u)/(n*a^2*factor)*N;

    statedot(4) = r*sin(u)/(n*a^2*factor*sin(i))*N;

    statedot(5) = factor/(n*a*e)*(-cos(nu)*R + sin(nu)*(1+r/p)*T)...
        -r*cot(i)*sin(u)/h*N;

    statedot(6) = h/r^2 + 1/(e*h)*p*cos(nu)*R-(p+r)*sin(nu)*T;

end