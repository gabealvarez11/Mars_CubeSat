function [value, isterminal, direction] = reachedMars(t,oe,mu,aMars)
    rv = equi2rv(oe,mu);
    r = norm(rv(1:3));
    
    value = (r>=aMars);
    isterminal = 1;
    direction = 0;
end