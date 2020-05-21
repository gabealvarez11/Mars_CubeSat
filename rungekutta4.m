function y = rungekutta4(ydot,y0,tspan)
    y = zeros(length(tspan),length(y0));
    y(1,:) = y0;
    for i = 2:length(tspan)
       h = tspan(i)-tspan(i-1);
       k1 = h*ydot(tspan(i-1),y(i-1,:).');
       k2 = h*ydot(tspan(i-1)+h/2,y(i-1,:).'+k1/2);
       k3 = h*ydot(tspan(i-1)+h/2,y(i-1,:).'+k2/2);
       k4 = h*ydot(tspan(i-1)+h,y(i-1,:).'+k3);
       y(i,:) = y(i-1,:).' + (1/6)*(k1+2*k2+2*k3+k4);
    end    
end