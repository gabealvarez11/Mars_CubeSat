function E = EfromM(M,e,ep)
    if e>= 0.5
        E = pi;
    else
        E = M;
    end
    d = delta(E,e,M);
    n = 0;
    while abs(d) > ep
       E = E+d;
       d = delta(E,e,M);
       n = n + 1;
       if n == 20
           break;
       end
    end
    
    function delta = delta(Ei,e,M)
        delta = (Ei - e*sin(Ei) - M)/(e*cos(Ei) - 1);
    end
end