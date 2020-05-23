function inertial_rv = kep2rv(state, mu)
% KEP2RV  Converts Keplerian elements to intertial position, velocity.
%   
%   INPUT: state = [a e i raan arg_periapsis nu], in [km, rad]
%          mu [km^3/s^2]
%   OUTPUT: inertial_rv = [r v]', in [km, km/s]

    a = state(1);
    e = state(2);
    i = rad2deg(state(3));
    raan = rad2deg(state(4));
    arg_periapsis = rad2deg(state(5));
    nu = rad2deg(state(6));
    
    % Semi-latus rectum.
    p = a * (1 - e^2);
    
    % Distance from center of the central body.
    r = p / (1 + e * cosd(nu));
    
    % Construct the position vector in the PQW frame.
    r_P = r * cosd(nu);
    r_Q = r * sind(nu);
    r_W = 0;
    
    r_PQW = [r_P; r_Q; r_W];
    
    % Find the speed.
    specific_energy = - mu / (2 * a);
    v = sqrt(2 * (specific_energy + mu / r));
    
    % Construct the velocity vector in the PQW frame.
    h = sqrt(mu * p);
    h_PQW = [0; 0; h];
        
    theta = asind(h / (r * v));
    
    v_PQW = (r * v * cosd(theta) * r_PQW + cross(h_PQW, r_PQW)) / r ^ 2;
    
    if i == 0
        if e == 0 % Equatorial and circular orbit: raan, arg_periapsis undefined.
            PQW_R_IJK = eye(3);
        else % Equatorial but non-circular: 
            PQW_R_IJK = rotz(arg_periapsis);
        end
    elseif e == 0 % Circular but non-equatorial orbit: arg_periapsis undefined.
        PQW_R_IJK = rotz(raan)* rotx(i);
    else % Non-circular, non-equatorial: all oe defined.
        PQW_R_IJK = rotz(raan) * rotx(i) * rotz(arg_periapsis);
    end
    
    r_IJK = PQW_R_IJK * r_PQW;
    v_IJK = PQW_R_IJK * v_PQW;
    
   inertial_rv = [r_IJK; v_IJK];
end

% function [r_eci,v_eci] = kep2ECI(oe, mu)
%     a = oe(1); e = oe(2); i = oe(3); raan = oe(4); w = oe(5); M = oe(6);
%     if e == 0
%        w = 0; 
%     end
%     if i == 0
%        raan = 0; 
%     end
%     %Calculate additional useful orbit parameters
%     nu = nufromE(EfromM(M,e,1e-10),e);
%     p = a*(1-e^2);
%     r = p/(1+e*cos(nu));
%     v = sqrt((2*mu/r) - (mu/a));
%     h = sqrt(mu * p);
%     
%     %Calculate angle between r and v vector
%     gamma = asin(h/(r*v));
%     %Calculate position and velocity vectors
%     rotm = rotz(rad2deg(raan))*rotx(rad2deg(i))*rotz(rad2deg(w));
%         % ... Rotation matrix from PQW to IJK 
%     r_vec_kep = [r*cos(nu); r*sin(nu); 0]; % PQW coordinates
%     v_vec_kep = (1/r^2).*((v*r*cos(gamma)).*r_vec_kep - cross(r_vec_kep, [0;0;h])); %calculate v in PQW from h,r,norm(v) and gamma
%     r_eci = (rotm*(r_vec_kep));
%     v_eci = (rotm*(v_vec_kep));
% end