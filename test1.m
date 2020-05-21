clearvars; clc; close all; set(0,'DefaultFigureWindowStyle','docked');
year2sec= 31556952;
Isp = 1200*9.81e-3; %decent EP, km/s
F = 16e-6*[0;1;0]; %kg*km/s2 = N/1000
Fs = (1:10)*1e-6;
F_dt = year2sec*(0.135); %sec
m_struct = 15; %kg
aEarth = 149.60e6; %km
mu_sun = 1.32712440018e11; %km3/s2
n = sqrt(mu_sun/aEarth^3); %rad/sec
mu_mars = 4.282837e4; %km3/s2

oe0 = [a; 0; 0; 0; 0; 0]; %start at Earth circular orbit (parabolic departure from Earth)

[dt,v_arr] = run_trial(oe0,F,F_dt,m_struct,Isp); dt = dt/year2sec;

% v_arr2 = zeros(length(F_dt),length(Fs));
% dt2 = zeros(length(F_dt),length(Fs));
% for j = 1:length(Fs)
%     [dt,v_arr] = run_trial(oe0,Fs(j)*[0;1;0],F_dt,m_struct,Isp);
%     v_arr2(:,j) = v_arr;
%     dt2(:,j) = dt;
% end
% dt2 = dt2/86400;
% F_dt = F_dt/86400;
% [DT,F] = meshgrid(F_dt,Fs);
% %% Plotting
% figure(); hold on;
% contourf(DT,F*1e6,dt2.');
% xlabel('Time of Powered Flight, days'); ylabel('Thrust, uN');
% title('Time-of-Flight to Mars, days'); colorbar;
% set(gca,'color','k');
% figure(); hold on; 
% contourf(DT,F*1e6,v_arr2.');
% xlabel('Time of Powered Flight, days'); ylabel('Thrust, uN');
% title('Arrival Velocity Relative to Mars, km/s'); colorbar;
% set(gca,'color','k');