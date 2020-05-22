 function PlotMe(hObject,event,datas,dts,v_arrs,prop_names,errs,mu_sun)
    mousePos =event.IntersectionPoint;
    F_dt = hObject.XData;
    Y = hObject.YData;
	[~,i] = min(abs(mousePos(1)-F_dt(1,:)));
    [~,j] = min(abs(mousePos(2)-Y(:,1)));
    n = get(gcf,'Number');
    p = floor((n-1)/3)+1;
    oe_out = datas{p,j};
    oe_out = oe_out{i};
    r_eci = zeros(size(oe_out,1),3);
    for k = 1:length(oe_out)
       curr_rv = equi2rv(oe_out(k,:),mu_sun);
       r_eci(k,:) = curr_rv(1:3);
    end
    ax = gca;
    h = figure('Name',[prop_names{p} ' ' 'Trajectory']); hold on;
    set(h,'WindowStyle','normal');
    title({strcat(ax.XLabel.String, ": ", num2str(F_dt(1,i)), " | ", ax.YLabel.String, ": ", num2str(Y(j,1))),...
       strcat("TOF, days: ", num2str(dts(p,i,j))," | V_{arr}, km/s: ",num2str(v_arrs(p,i,j))," | Error: ",num2str(errs(p,i,j)))});
    xlabel('Heliocentric X (km)'); ylabel('Heliocentric Y (km)');
    aMars = 227939200; %km
    aEarth = 149.60e6; %km
    viscircles([0 0],aEarth,'Color','b','LineStyle','--','LineWidth',1);
    viscircles([0 0],aMars,'Color','r','LineStyle','--','LineWidth',1);
    tstep = dts(p,i,j)/size(r_eci,1);
    if errs(p,i,j) == 3 || errs(p,i,j) == 1
        plot(r_eci(:,1),r_eci(:,2),'g','DisplayName','Continuous Thrust');
    elseif errs(p,i,j) ~= 2
        cuti = floor(F_dt(1,i)/tstep);
        plot(r_eci(1:cuti,1),r_eci(1:cuti,2),'g','DisplayName','Continuous Thrust');
        plot(r_eci(cuti:end,1),r_eci(cuti:end,2),'k','DisplayName','Ballistic');
    else
        text(0,0,"Simulation Timeout ERROR");
        plot(r_eci(:,1),r_eci(:,2),'g','DisplayName','Trajectory');
    end
    legend; axis equal;
 end
