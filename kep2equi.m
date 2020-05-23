function equi = kep2equi(kep)
% KEP2EQUI  Converts Keplerian elements to equinoctial elements.
%   
%   INPUT: kep = [a e i raan arg_periapsis nu], in [rad]
%   OUTPUT: equi = [p f g h k L], in [rad]
    
    a = kep(1);
    e = kep(2);
    i = kep(3);
    raan = kep(4);
    arg_periapsis = kep(5);
    nu = kep(6);
    
    equi = zeros(6,1);
    
    equi(1) = a * (1-e^2);
    equi(2) = e * cos(arg_periapsis + raan);
    equi(3) = e * sin(arg_periapsis + raan);
    equi(4) = tan(i/2) * cos(raan);
    equi(5) = tan(i/2) * sin(raan);
    equi(6) = arg_periapsis + raan + nu;
end