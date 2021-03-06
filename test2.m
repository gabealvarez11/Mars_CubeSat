clearvars; clc; close all; set(0,'DefaultFigureWindowStyle','docked');
year2sec= 31556952;
prop_names = {'FEEP', 'Hall Thruster', 'Electrospray'};
Isps = [3000; 1200; 800]*9.81e-3; %decent EP, km/s
Fs = [.35; 16; 0.03]*1e-6; %kg*km/s2 = N/1000
Fangle = deg2rad(0:5:90);
F_dt = year2sec*(0.25:0.25:9.75); %sec
m_struct = [1; 25; 1]; %kg
aEarth = 149.60e6; %km
mu_sun = 1.32712440018e11; %km3/s2
n = sqrt(mu_sun/aEarth^3); %rad/sec
mu_mars = 4.282837e4; %km3/s2

oe0 = [aEarth; 0; 0; 0; 0; 0]; %start at Earth circular orbit (parabolic departure from Earth)

v_arrs = zeros(length(prop_names),length(F_dt),length(Fangle));
dts = zeros(length(prop_names),length(F_dt),length(Fangle));
errs = zeros(length(prop_names),length(F_dt),length(Fangle));
datas = cell(length(prop_names),length(Fangle));
for k = 1:length(prop_names)
    for j = 1:length(Fangle)
        [dt,v_arr,errors,data] = run_trial(oe0,Fs(k)*[sin(Fangle(j));cos(Fangle(j));0],F_dt,m_struct(k),Isps(k));
        v_arrs(k,:,j) = v_arr;
        dts(k,:,j) = dt;
        errs(k,:,j) = errors;
        datas{k,j} = data;
    end
end
dts = dts/86400;
F_dt = F_dt/86400;
%% Plotting
for i = 1:length(prop_names)
    [DT,FA] = meshgrid(F_dt,Fangle);
    dtp = squeeze(dts(i,:,:));
    v_arrp = squeeze(v_arrs(i,:,:));
    errp = squeeze(errs(i,:,:));
    figure('Name',[prop_names{i} ' ' 'TOF']); hold on;
    h1 = surf(DT,rad2deg(FA),dtp.');
    xlabel('Time of Powered Flight, days'); ylabel('Angle, deg');
    title('Time-of-Flight to Mars, days'); colorbar;
    set(gca,'color','k');
    figure('Name',[prop_names{i} ' ' 'V Arr']); hold on; 
    h2 = surf(DT,rad2deg(FA),v_arrp.');
    xlabel('Time of Powered Flight, days'); ylabel('Angle, deg');
    title('Arrival Velocity Relative to Mars, km/s'); colorbar;
    set(gca,'color','k');
    figure('Name',[prop_names{i} ' ' 'Errors']); hold on; 
    h3 = surf(DT,rad2deg(FA),errp.');
    xlabel('Time of Powered Flight, days'); ylabel('Angle, deg');
    title('Error Codes'); colorbar;
    set(gca,'color','k');
    h1.ButtonDownFcn = @(hO,e) PlotMe(hO,e,datas,dts,v_arrs,prop_names,errs,mu_sun);
    h2.ButtonDownFcn = @(hO,e) PlotMe(hO,e,datas,dts,v_arrs,prop_names,errs,mu_sun);
    h3.ButtonDownFcn = @(hO,e) PlotMe(hO,e,datas,dts,v_arrs,prop_names,errs,mu_sun);
end

