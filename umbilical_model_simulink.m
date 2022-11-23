function [r,v,a] = umbilical_model_simulink(t,cable_length,v_ship,current,waves,Ts,r_i,r_ip1,r_im1,v_i)
    

    % Define the parameter values used in the model
    
    n = 20; % # of segments
    k = n+1; % # of nodes
    g = [0,0,9.81]; % m/s^2, gravitational acceleration
    rho = 1020; % kg/m^3, density of seawater
    d = 0.01; % m, diameter of umbilical
    vc = current; % m/s, water current
    E = 7.2*10^10; % Pa, kg/(m*s^2), Young's modulus
    lc = cable_length; % m, umbilical cable length
    m_km = 320; % kg/km, weight per km umbilical cable
    mc = m_km * lc/1000; % kg, total weight for umbilical cable
    Ct = 0.2; % Tangential drag coefficient
    Cn = 0.4; % Normal drag coefficient
    l0 = lc/n; % initial length of all segments
    m = eye(3)*mc/n; % kg, weight per segment
    
    % Initializing n by 3 matrices for position, velocity and acceleration
    
    r = zeros(n,3);
    v = zeros(n,3);
    a = zeros(n,3);
    
    % Set the velocity of the first node to be equal to the ships velocity
    
    v(1,:) = v_ship;
    
    % Set the position of the first node to be equal to the ships estimated
    % position.
    % ATTENTION this a forward Euler prediction and it is only
    % used if the position of the ship is unavailable. Otherwise add an
    % input in the function and set up the first node position to be equal
    % to the ships position.
    
    r(1,:) = r(1,:)+v(1,:)*Ts;


    %% Runge-Kutta fourth order algorithm for all nodes, except node 1 and n+1
    for j=2:n
        
        % Set current node position equal to the prediction from last step
        
        r(j,:) = r_ip1(j,:);
        
        % Use the system of first order ODEs, from the Runge-Kutta method,
        % in order to correct the prediction and calculate the velocity and
        % acceleration.
        
        K1r = Ts*f1(t,r_i(j,:),v_i(j,:));
        K1v = Ts*f2(t,r_i(j,:),v_i(j,:),r_ip1(j,:),r_im1(j,:),r_i(j-1,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K2r = Ts*f1(t+Ts/2,r_i(j,:)+K1r/2,v_i(j,:)+K1v/2);
        K2v = Ts*f2(t+Ts/2,r_i(j,:)+K1r/2,v_i(j,:)+K1v/2,r_ip1(j,:),r_im1(j,:),r_i(j-1,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K3r = Ts*f1(t+Ts/2,r_i(j,:)+K2r/2,v_i(j,:)+K2v/2);
        K3v = Ts*f2(t+Ts/2,r_i(j,:)+K2r/2,v_i(j,:)+K2v/2,r_ip1(j,:),r_im1(j,:),r_i(j-1,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K4r = Ts*f1(t+Ts,r_i(j,:)+K3r/2,v_i(j,:)+K3v/2);
        K4v = Ts*f2(t+Ts,r_i(j,:)+K3r/2,v_i(j,:)+K3v/2,r_ip1(j,:),r_im1(j,:),r_i(j-1,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        
        % Updating the r(j,:) and v(j,:) values using the Runge-Kutte
        % method.
        
        r(j,:) = r_i(j,:) + 1/6*(K1r+2*K2r+2*K3r+K4r);
        v(j,:) = v_i(j,:) + 1/6*(K1v+2*K2v+2*K3v+K4v);
        a(j,:) = K1v;
        
        % Calculate the sum of forces.
        
        if norm(r(j-1,:)-r(j,:)) > l0 && all(abs(buoyancy(d,g,l0,rho)-gravity(g,m)-drag(v_i(j,:),r_ip1(j,:),r_i(j,:),r_im1(j,:),Cn,Ct,d,rho,vc)) > abs(tension(r_i(j,:),r_i(j-1,:),d,E,l0)-tension(r_i(j,:),r_i(j-1,:),d,E,l0)))
            
            % Do not apply spherical boundary conditions, since the sum of
            % the forces and the tension will be applied. 
            % This only applies if the segment is stretched and the segment
            % is supposed to stretch.
            
        elseif norm(r(j-1,:)-r(j,:)) < l0 || norm(r(j-1,:)-r(j,:)) > l0
            
            % Apply the spherical boundary conditions, since the segment is
            % not supposed to stretch.
            
            P = r(j,:) - r(j-1,:);
            Q = (l0/norm(P))*P;
            r(j,:) = Q + r(j-1,:);
        end
    end
end
%% Functions for equation of motion and functions for Runge-Kutta method

% Function for the first, 1st order ODE for the RK method

function f1 = f1(t,r_i,v_i)
    f1 = v_i;
end

% Function for the second, 1st order ODE for the RK method

function f2 = f2(t,r_i,v_i,r_ip1,r_im1,r_jm1,g,m,d,E,l0,Cn,Ct,rho,vc,waves)
    if abs(buoyancy(d,g,l0,rho)-gravity(g,m)-drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc)) > abs(tension(r_i,r_jm1,d,E,l0)-tension(r_i,r_jm1,d,E,l0))
        f2 = (buoyancy(d,g,l0,rho)-gravity(g,m)-drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc)+waves+tension(r_i,r_jm1,d,E,l0)-tension(r_i,r_jm1,d,E,l0))/(mass(d,l0,m,rho));
    else
        f2 = (buoyancy(d,g,l0,rho)-gravity(g,m)-drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc)+waves)/(mass(d,l0,m,rho));
    end
end

% Function for drag force, using the tangential and normal components

function Fd = drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc)
    t_i = (r_ip1-r_im1)/norm(r_ip1-r_im1);
    D_t = rho*pi*d/4*Ct*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc').*t_i).*(v_i-vc').*t_i;
    D_n = rho*pi*d/4*Cn*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc')-(v_i-vc').*t_i)*((v_i-vc')-(v_i-vc').*t_i);
    Fd = D_t + D_n;
end 

% Function for tension force, currently not used

function Ft = tension(r_i,r_jm1,d,E,l0)
    Ft = pi*d^2*E*(norm(r_jm1-r_i)-l0)*(r_jm1-r_i)/(4*l0*norm(r_jm1-r_i));
    
    % If its compression that occurs, then the tension should be 0
    
    if norm(r_i-r_jm1) < l0
        Ft = [0,0,0];
    end
end

% Function for gravitational force

function Fg = gravity(g,m)
    Fg = g*m;
end

% Function for buoyancy force

function Fb = buoyancy(d,g,l0,rho)
    Fb = g*rho*l0*d^2/4;
end

% Function for total mass of node

function Fm = mass(d,l0,m,rho)
    Fm = m+eye(3)*rho*pi*d^2/4*l0;
end
