clearvars; clc; close all; set(0,'DefaultFigureWindowStyle','docked');
year2sec= 31556952;
prop_names = {'FEEP', 'Hall Thruster', 'Electrospray'};
Isps = [3000; 1200; 800]*9.81e-3; %decent EP, km/s
Fs = [.35 16 0.03]*1e-6*[0;1;0]; %kg*km/s2 = N/1000
F_dt = year2sec*(0.5:0.5:5); %sec
m_struct = [15; 15; 15]; %kg
aEarth = 149.60e6; %km
mu_sun = 1.32712440018e11; %km3/s2
n = sqrt(mu_sun/aEarth^3); %rad/sec
mu_mars = 4.282837e4; %km3/s2

oe0 = [a; 0; 0; 0; 0; 0]; %start at Earth circular orbit (parabolic departure from Earth)

[dt,v_arr] = run_trial(oe0,F,F_dt,m_struct,Isp); dt = dt/year2sec;

v_arrs = zeros(length(prop_names),length(F_dt),length(Fs));
dts = zeros(length(prop_names),length(F_dt),length(Fs));
errs = zeros(length(prop_names),length(F_dt),length(Fs));
for k = 1:length(prop_names)
    for j = 1:length(Fs)
        [dt,v_arr,errors] = run_trial(oe0,Fs(j)*[0;1;0],F_dt,m_struct,Isp);
        v_arrs(k,:,j) = v_arr;
        dts(k,:,j) = dt;
        errs(k,:,j) = errors;
    end
end
dts = dts/86400;
F_dt = F_dt/86400;
%% Plotting
for i = 1:length(prop_names)
    [DT,F] = meshgrid(F_dt(i,:,:),Fs(i,:,:));
    figure('Name',[prop_names{i} ' ' 'TOF']); hold on;
    contourf(DT,F*1e6,dts(i,:,:).');
    xlabel('Time of Powered Flight, days'); ylabel('Thrust, mN');
    title('Time-of-Flight to Mars, days'); colorbar;
    set(gca,'color','k');
    figure('Name',[prop_names{i} ' ' 'TOF']); hold on; 
    contourf(DT,F*1e6,v_arrs(i,:,:).');
    xlabel('Time of Powered Flight, days'); ylabel('Thrust, uN');
    title('Arrival Velocity Relative to Mars, km/s'); colorbar;
    set(gca,'color','k');
end