function kep = equi2kep(equi)
% KEP2EQUI  Converts Keplerian elements to equinoctial elements.
%   
%   INPUT: equi = [p f g h k L], in [rad]
%   OUTPUT: kep = [a e i raan arg_periapsis nu], in [rad]

    p = equi(1); f = equi(2); g = equi(3);
    h = equi(4); k = equi(5); L = equi(6);
    
    kep = zeros(6,1);
    
    kep(1) = p / (1 - f^2 - g^2);
    kep(2) = sqrt(f^2 + g^2);
    kep(3) = atan2(2*sqrt(h^2+k^2),1-h^2-k^2);
    kep(4) = atan2(k,h);
    kep(5) = atan2(g*h-f*k,f*h+g*k);
    kep(6) = L - atan(g/f);
end
    