function nu = nufromE(E,e)
    E = mod(E,2*pi);
    nu = acos((cos(E)-e)/(1-e*cos(E)));
    if ((E>=0 && E<=pi) && ~(nu>=0 && nu<=pi)) || ((E>pi && E<2*pi) && ~(nu>pi && nu<2*pi))
        nu = 2*pi - nu;
    end
end