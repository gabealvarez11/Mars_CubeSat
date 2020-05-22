function [dt, v_arr,errors] = run_trial(oe0,F,F_dt,m_struct,Isp)
% Propagator for along-track constant thrust, using orbital elements and GVEs to provide state change
% INPUTS:
% ---------------------------------------
% tstep: size of time steps in propogation, in seconds
% oe0: [a, e, i, raan, w, M0], all angles in radians
% f_dir: unit vector with components [fR fT fN] (radial tangential normal) 
%       giving the direction of the continuous thrust vector
% f_dt: duration of thrust period, seconds
% m_struct: structural mass (dry mass) of spacecraft, in kg
% Isp: specific impulse of propulsion system, in m/s (Isp*g0)
% opts: length-4 array whose non-zero elements reset the defaults 
%       options [mu_sun,tmax,Mars_escape_vel,close_2_Mars_tol]
%       mu_sun: center body grav parameter, km3/s2
%       tmax: cut-off time for simulation, sec
%       Mars_escape_vel: max relative velocity to Mars to be captured, km/s
%       close_2_Mars_tol: dist from Mars at which to cut off simulation, km
%
% OUTPUTS:
% ---------------------------------------
% dt: time-of-flight from Earth departure to Mars capture arrival, in sec
% v_ex: velocity relative to Mars upon arrival, in km/s
% error: string explaining why simulation cut off ('tooFast' or 'tooLong',
%           indicating too much excess velocity at arrival or a failure to
%           reach Mars orbital radius before tmax was reached. Empty
%           if trajectory was successful).
    %% Constants/Defaults
    year2sec= 31556952;
    aMars = 227939200; %km
    mu_sun = 1.32712440018e11; %km3/s2
    tmax = max(10*year2sec,max(F_dt)); %sec, max duration of simulation before cut-off
    tstep = tmax/1e5;
    %% Calculate Values
    vMars = sqrt(mu_sun/aMars);
    mdot = norm(F)/Isp;
    %% Initialize Variables
    sim_opts = odeset('AbsTol',1e-6,'RelTol',1e-9,'Events',@(t,oe) reachedMars(t,oe,mu_sun,aMars)); % set propogator tolerance
    dt = zeros(size(F_dt)); %initialize TOF vector
    v_arr = zeros(size(F_dt)); %initialize excess velocity vector
    errors = zeros(size(F_dt));
    
    E0 = EfromM(oe0(6),oe0(2),1E-10);
    oe0(6) = nufromE(E0,oe0(2));
    oe0 = kep2equi(oe0);
    
    %% Iterate Thrust Durations
    for j  = 1:length(F_dt)
        t_cutoff = F_dt(j); % get t_cutoff for this run
        m0 = m_struct+mdot*t_cutoff; %get starting wet mass for this run
        sim_error = 0; %initialize simulation error variable
        %[~, oe_lowthrust,tof,~,~] = ode113(@getdot_lowthrust,0:tstep:t_cutoff,oe0,sim_opts1);
        [~, oe_lowthrust,tof,~,~] = ode113(@(t,state) GVEs_equi(t,state,mu_sun,F,m0,mdot),0:tstep:t_cutoff,oe0,sim_opts);

        curr_e = norm(oe_lowthrust(end,2:3));
        rapogee = oe_lowthrust(end,1)/(1-curr_e);
        
        if isempty(tof) && rapogee >= aMars
            %Mball = oe_lowthrust(end,6);
            [~,oe_ballistic,tof,~,~] = ode113(@(t,state) getdot_ballistic(t,state,mu_sun),t_cutoff:tstep:tmax,oe_lowthrust(end,:).',sim_opts);
            oe_out = [oe_lowthrust;oe_ballistic];
        else
            if rapogee < aMars
                tof = -1;
                sim_error = 1; %won't reach Mars
            else
                sim_error = 3; % already there.
            end
            oe_out = oe_lowthrust;
        end
        if isempty(tof) && sim_error==0 
            tof = -1;
            sim_error= 2; %simulation timed out before event occured
        end
        
        curr_rv = equi2rv(oe_out(length(oe_out),:),mu_sun);
        vf = curr_rv(4:6);
        dt(j) = tof(1);

        Lf = oe_out(length(oe_out),6);
        
        v_arr(j) = norm(vf - vMars*[-sin(Lf); cos(Lf); 0]); %get relative speed to Mars (circular)
        errors(j) = sim_error;
        % Plot test
%             figure(); hold on;
%             r_eci = zeros(size(oe_out,1),3);
%             for k = 1:length(oe_out)
%                rv1 = equi2rv(oe_out(k,:),mu_sun);
%                r_eci(k,:) = rv1(1:3);
%             end
%             aEarth = 149.60e6; %km
%             viscircles([0 0],aEarth,'Color','b','LineStyle','--','LineWidth',1);
%             viscircles([0 0],aMars,'Color','r','LineStyle','--','LineWidth',1);
%             cuti = floor(t_cutoff/tstep);
%             if cuti>floor(tof/tstep) && sim_error ~= 2
%                 plot(r_eci(:,1),r_eci(:,2),'g','DisplayName','Continuous Thrust');
%             else
%                 plot(r_eci(1:cuti,1),r_eci(1:cuti,2),'g','DisplayName','Continuous Thrust');
%                 plot(r_eci(cuti+1:end,1),r_eci(cuti+1:end,2),'k','DisplayName','Ballistic');
%             end
%             axis equal; xlabel('Heliocentric X, km'); ylabel('Heliocentric Y, km'); 
%             title('Low-Thrust Earth to Mars Trajectory'); legend;
    end
    %% Filter Arrival Speeds for Errors, Excessive Velocity
    dt(errors~=0 & errors ~= 3) = NaN;
    v_arr(errors~=0 & errors ~= 3) = NaN;



    
    
end