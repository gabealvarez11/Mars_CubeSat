function [dt, v_arr] = run_trial(oe0,F,F_dt,m_struct,Isp)
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
    sim_opts1 = odeset('AbsTol',1e-6,'RelTol',1e-9,'Events',@reachedMars); % set propogator tolerance
    sim_opts2 = odeset('AbsTol',1e-6,'RelTol',1e-9,'Events',@reachedMars_noMulti);
    dt = zeros(size(F_dt)); %initialize TOF vector
    v_arr = zeros(size(F_dt)); %initialize excess velocity vector
    errors = zeros(size(F_dt));
    %% Iterate Thrust Durations
    for j  = 1:length(F_dt)
        t_cutoff = F_dt(j); % get t_cutoff for this run
        m0 = m_struct+mdot*t_cutoff; %get starting wet mass for this run
        sim_error = 0; %initialize simulation error variable
        [~, oe_lowthrust,tof,~,~] = ode113(@getdot_lowthrust,0:tstep:t_cutoff,oe0,sim_opts1);
        rapogee = oe_lowthrust(end,1)*(1+oe_lowthrust(end,2));
        if isempty(tof) && rapogee >= aMars
            Mball = oe_lowthrust(end,6);
            [~,oe_ballistic,tof,~,~] = ode113(@getdot_ballistic,t_cutoff:tstep:tmax,oe_lowthrust(end,:).',sim_opts2);
            oe_out = [oe_lowthrust;oe_ballistic];
        else
            if rapogee < aMars
                tof = -1;
                sim_error = 1; %won't reach Mars
            end
            oe_out = oe_lowthrust;
        end
        if isempty(tof) && sim_error==0 
            tof = -1;
            sim_error= 2; %simulation timed out before event occured
        end
        [~,vf] = kep2ECI(oe_out(length(oe_out),:),mu_sun);
        dt(j) = tof(1);
        Mf = oe_out(length(oe_out),6); %final mean anomaly at Mars capture
        ef = oe_out(length(oe_out),2);
        nuf = nufromE(EfromM(Mf,ef,1e-10),ef);
        uf = nuf + oe_out(length(oe_out),5);
        v_arr(j) = norm(vf - vMars*[-sin(uf); cos(uf); 0]); %get relative speed to Mars (circular)
        errors(j) = sim_error;
        % Plot test
            figure(); hold on;
            r_eci = zeros(size(oe_out,1),3);
            for k = 1:length(oe_out)
               [r1, ~] = kep2ECI(oe_out(k,:),mu_sun);
               r_eci(k,:) = r1;
            end
            aEarth = 149.60e6; %km
            viscircles([0 0],aEarth,'Color','b','LineStyle','--','LineWidth',1);
            viscircles([0 0],aMars,'Color','r','LineStyle','--','LineWidth',1);
            cuti = floor(t_cutoff/tstep);
            if cuti>floor(tof/tstep) && sim_error ~= 2
                plot(r_eci(:,1),r_eci(:,2),'g','DisplayName','Continuous Thrust');
            else
                plot(r_eci(1:cuti,1),r_eci(1:cuti,2),'g','DisplayName','Continuous Thrust');
                plot(r_eci(cuti+1:end,1),r_eci(cuti+1:end,2),'k','DisplayName','Ballistic');
            end
            axis equal; xlabel('Heliocentric X, km'); ylabel('Heliocentric Y, km'); 
            title('Low-Thrust Earth to Mars Trajectory'); legend;
    end
    errors
    %% Filter Arrival Speeds for Errors, Excessive Velocity
    dt(errors~=0) = NaN;
    v_arr(errors~=0) = NaN;
%% Helper Functions
    % Need f_R, f_T, f_N, mu_sun, r_mars.
    % oe = [a e i raan arg_periapsis M] in [deg].
    function oedot = getdot_lowthrust(t,oe)
        oedot = zeros(6,1);
        a = oe(1);
        e = oe(2);
        i = oe(3);
        raan = oe(4);
        arg_periapsis = oe(5);
        M = oe(6);

        % Convert to true anomaly.
        E = EfromM(M,e,1E-10);
        nu = nufromE(E,e);
        % Calculate useful orbital information.
        n = sqrt(mu_sun/a^3);
        beta = sqrt(1-e^2);
        p = a*(1-e^2);
        r = p/(1+e*cos(nu));
        
        %Calculate current mass, specific force
        m = m0 - mdot*t; %kg
        f_R = F(1)/m; % m/s2
        f_T = F(2)/m;
        f_N = F(3)/m;
        % Use the GVEs to model continuous thrust equations of motion.
        oedot(1) = 2/n/beta*(e*sin(nu)*f_R + p/r*f_T);

        oedot(2) = beta/n/a*(sin(nu)*f_R + ((e+cos(nu))/(1+e*cos(nu)))*f_T);

        oedot(3) = r*cos(arg_periapsis+nu)/(n*a^2*beta)*f_N;

        oedot(4) = 0; %r*sin(arg_periapsis+nu)/n/(a^2)/beta/sin(i)*f_n

        oedot(5) =beta/(n*a*e)*sin(nu)*(1+r/p)*f_T -beta*cos(nu)/(n*a*e)*f_R; %...
            %-r*cot(i)*sin(arg_periapsis+nu)/(n*a^2*beta)*f_N;

        oedot(6) = n + 1/(n*a^2*e)*(p*cos(nu)-2*e*r)*f_R - 1/(n*a^2*e)*(p+r)*sin(nu)*f_T;
    end

    function oedot = getdot_ballistic(t,oe)
        oedot = zeros(6,1);
        n = sqrt(mu_sun/oe(1)^3);
        oedot(6) = n;
        % Propagate nu. All other oe are constant.
    end
    
    function [value, isterminal, direction] = reachedMars(t,oe)
        a = oe(1);
        e = oe(2);
        M = oe(6);
        % Convert to true anomaly.
        E = EfromM(M,e,1E-10);
        nu = nufromE(E,e);
        % Calculate useful orbital information.
        p = a*(1-e^2);
        r = p/(1+e*cos(nu));
       value = (r>=aMars);
       isterminal = 1;
       direction = 0;
    end

    function [value, isterminal, direction] = reachedMars_noMulti(t,oe)
        a = oe(1);
        e = oe(2);
        M = oe(6);
        % Convert to true anomaly.
        E = EfromM(M,e,1E-10);
        nu = nufromE(E,e);
        % Calculate useful orbital information.
        p = a*(1-e^2);
        r = p/(1+e*cos(nu));
       value = [(r>=aMars); (abs(M-Mball)>2*pi)];
       isterminal = [1; 1];
       direction = [0; 0];
       if value(2)
          sim_error = 2; 
       end
    end
end