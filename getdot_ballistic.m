function oedot = getdot_ballistic(t,oe,mu)
    p = oe(1); f = oe(2); g = oe(3);
    L = oe(6);
    
    oedot = zeros(6,1);
    w = 1 + f * cos(L) + g * sin(L);
    
    oedot(6) = sqrt(mu*p)*(w/p)^2;
end