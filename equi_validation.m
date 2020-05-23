clc; clear all; close all; set(0,'DefaultFigureWindowStyle','docked');

Isp = 1200*9.81e-3; % [km/s]
thrust = 16*1e-8; % [kg*km/s^2] = [mN]

%thrust = 9.2e-6; % override for testing

mu_E = 3.986e5; % [km^3/s^2]

F = zeros(4,3);
F(1,1) = thrust;
F(2,2) = thrust;
F(3,3) = thrust;
dir = ["Radial Thrust" "Tangential Thrust" "Normal Thrust" "No Perturbation"];
dir = convertStringsToChars(dir);

%F = [0 thrust 0];
%dir = {'Tangential Thrust'};

m0 = 10; % [kg]
mdot = norm(thrust)/Isp; % [kg/s]

options = odeset('AbsTol',1e-6,'RelTol',1e-9);

a = 7000; % [km]
e = 0.2;
i = deg2rad(30); % all angles in [rad]
raan = deg2rad(40);
arg_periapsis = deg2rad(50);
nu = deg2rad(160);

oe0_kep = [a e i raan arg_periapsis nu];
oe0_equi = kep2equi(oe0_kep);

oe_kep_check = equi2kep(oe0_equi);
%oe0_kep(6) = nu2M(e,nu);

rv0_kep = kep2rv(oe0_kep,mu_E);
rv0_equi = equi2rv(oe0_equi,mu_E);

duration = 30 * 86400; % [s]
n = 1E4;
tstep = duration / n;
n = n + 1;

numCurves = length(F(:,1));

figure(1); hold on;
%figure(2); hold on;

for j=1:numCurves
    
    [t_out, kep_vallado_out] = ode113(@(t,state) GVEs_kep_vallado(t,state,mu_E,F(j,:),m0,mdot),0:tstep:duration,oe0_kep,options);
    [~, kep_esoc_out] = ode113(@(t,state) GVEs_kep_esoc(t,state,mu_E,F(j,:),m0,mdot),0:tstep:duration,oe0_kep,options);

    [~, equi_out] = ode113(@(t,state) GVEs_equi(t,state,mu_E,F(j,:),m0,mdot),0:tstep:duration,oe0_equi,options);
    
    equi_rv = zeros(n,6);
    kep_vallado_rv = zeros(n,6);
    kep_esoc_rv = zeros(n,6);
    
    deviation = zeros(n,1);
    fractional_error = zeros(n,1);
    
    % [EQUI KEP]
    specific_energy = zeros(n,2);
    
    for k=1:n
        equi_rv(k,:) = equi2rv(equi_out(k,:),mu_E);
        
        kep_vallado_rv(k,:) = kep2rv(kep_vallado_out(k,:),mu_E);
        kep_esoc_rv(k,:) = kep2rv(kep_esoc_out(k,:),mu_E);
        
        % Use Vallado as benchmark.
        %deviation(k) = norm(equi_rv(k,1:3) - kep_vallado_rv(k,1:3));

        % Or, use ESOC as benchmark.
        deviation(k) = norm(equi_rv(k,1:3) - kep_esoc_rv(k,1:3));
        
        % Compare the two Keplerian methods.
        %deviation(k) = norm(kep_vallado_rv(k,1:3) - kep_esoc_rv(k,1:3));
        
        %specific_energy(k,1) = norm(equi_rv(k,4:6))^2/2 - mu_E/norm(equi_rv(k,1:3));
        %specific_energy(k,2) = norm(kep_esoc_rv(k,4:6))^2/2 - mu_E/norm(kep_esoc_rv(k,1:3));
        
%         fractional_error(k) = norm(equi_rv(k,1:3) - kep_esoc_rv(k,1:3))/...
%             norm(equi_rv(k,1:3));
    end
    
    %figure('Name',[dir{j} ', mag'])
    %figure('Name',dir{j})
    
    figure(1);
    plot(t_out/86400,deviation,"DisplayName", dir{j});
    title("Accumulated Propagation Error, Vallado Keplerian Against ESOC Keplerian")
    legend;
    xlabel("t [days]");
    ylabel("Magnitude of Position Error [km]")
    
    % Check orbital energy.
    
%     figure(2);
%     plot(t_out/86400,specific_energy(:,1),"DisplayName", dir{j});
%     title("Orbit Raising Effect of Continuous Thrust")
%     legend;
%     xlabel("t [days]");
%     ylabel("Specific Energy [km^2/s^2]")
    
    % Fractional position error.
    
%     figure('Name',[dir{j} ', %'])
%     plot(t_out/86400,fractional_error,"DisplayName", dir{j});
%     legend;
%     xlabel("t [days]");
%     ylabel("Fractional Position Error, r_{error} / r")
    
    % Plotting slices of trajectories.
    
%     figure('Name',['Kep Traj. for ' dir{j} ', xy'])
%     plot(kep_rv(:,1),kep_rv(:,2),".");
%     xlabel("x"); ylabel("y"); axis equal;
%     
%     figure('Name',['Kep Traj. for ' dir{j} ', xz'])
%     plot(kep_rv(:,1),kep_rv(:,3),".");
%     xlabel("x"); ylabel("z"); axis equal;
%     
%     figure('Name',['Equi. Traj. for ' dir{j} ', xy'])
%     plot(equi_rv(:,1),equi_rv(:,2),".");
%     xlabel("x"); ylabel("y"); axis equal;  
%     
%     figure('Name',['Equi. Traj. for ' dir{j} ', xz'])
%     plot(equi_rv(:,1),equi_rv(:,3),".");
%     xlabel("x"); ylabel("z"); axis equal;
end

