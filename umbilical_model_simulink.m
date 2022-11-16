function [r,v,a] = umbilical_model(t,c_length,segments,v_ship,current,waves,Ts,r_i,r_ip1,r_im1,v_i)
    
    % Define the parameter values
    dt = Ts; % s, timestep
    n = 20; % # of segments
    k = n+1; % # of nodes
    g = [0,0,9.81]; % m/s^2, gravitational acceleration
    rho = 1020; % kg/m^3, density of seawater
    d = 0.01; % m, diameter of umbilical
    vc = current; % m/s, water current
    E = 7.2*10^10; % Pa, kg/(m*s^2), Young's modulus
    lc = c_length; % m, umbilical cable length
    m_km = 320; % kg/km, weight per km umbilical cable
    mc = m_km * lc/1000; % kg, total weight for umbilical cable
    Ct = 0.2; % Tangential drag coefficient
    Cn = 0.4; % Normal drag coefficient
    l0 = lc/n; % initial length of all segments
    m = eye(3)*mc/n; % kg, weight per segment
    
    % Initializing three n by 3 matrices
    r = zeros(n,3);
    v = zeros(n,3);
    a = zeros(n,3);
    v(1,:) = v_ship;

    %% Runge-Kutta fourth order algorithm
    h = dt;
    for j=2:n
        r(j,:) = r_ip1(j,:);
        K1r = h*f1(t,r_i(j,:),v_i(j,:));
        K1v = h*f2(t,r_i(j,:),v_i(j,:),r_ip1(j,:),r_im1(j,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K2r = h*f1(t+h/2,r_i(j,:)+K1r/2,v_i(j,:)+K1v/2);
        K2v = h*f2(t+h/2,r_i(j,:)+K1r/2,v_i(j,:)+K1v/2,r_ip1(j,:),r_im1(j,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K3r = h*f1(t+h/2,r_i(j,:)+K2r/2,v_i(j,:)+K2v/2);
        K3v = h*f2(t+h/2,r_i(j,:)+K2r/2,v_i(j,:)+K2v/2,r_ip1(j,:),r_im1(j,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves);
        K4r = h*f1(t+h,r_i(j,:)+K3r/2,v_i(j,:)+K3v/2);
        K4v = h*f2(t+h,r_i(j,:)+K3r/2,v_i(j,:)+K3v/2,r_ip1(j,:),r_im1(j,:),g,m,d,E,l0,Cn,Ct,rho,vc,waves');
        
        % Updating the r(j,:,i+1) value using the Runge-Kutta method
        r(j,:) = r_i(j,:) + 1/6*(K1r+2*K2r+2*K3r+K4r);
        v(j,:) = v_i(j,:) + 1/6*(K1v+2*K2v+2*K3v+K4v);
        a(j,:) = K1v;

        % Adjust for spherical boundary conditions
        if norm(r(j-1,:)-r(j,:)) < l0 || norm(r(j-1,:)-r(j,:)) > l0
            fa = fa_s(l0,r(j,1),r(j,2),r(j,3),r(j-1,1),r(j-1,2),r(j-1,3));
            r(j,:) = r(j,:)*fa;
        end
    end
end
%% Functions for equation of motion and functions for Runge-Kutta method

% Function for the first, 1st order ODE for the RK method
function f1 = f1(t,r_i,v_i)
    f1 = v_i;
end

% Function for the second, 1st order ODE for the RK method
function f2 = f2(t,r_i,v_i,r_ip1,r_im1,g,m,d,E,l0,Cn,Ct,rho,vc,waves)
    f2 = (buoyancy(d,g,l0,rho)-gravity(g,m)-drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc))/(mass(d,l0,m,rho));
end

% Function for drag force, using the tangential and normal components
function Fd = drag(v_i,r_ip1,r_i,r_im1,Cn,Ct,d,rho,vc)
    t_i = (r_ip1-r_im1)/norm(r_ip1-r_im1);
    D_t = rho*pi*d/4*Ct*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc').*t_i).*(v_i-vc').*t_i;
    D_n = rho*pi*d/4*Cn*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc')-(v_i-vc').*t_i)*((v_i-vc')-(v_i-vc').*t_i);
    Fd = D_t + D_n;
end 

% Function for tension force, currently not used
function Ft = tension(r_i,r_ip1,d,E,l0)
    Ft = pi*d^2*E*(norm(r_ip1-r_i)-l0)*(r_ip1-r_i)/(4*l0*norm(r_ip1-r_i));
    if norm(r_ip1-r_i) < l0
        Ft = 0; 
    end
end

% Function for gravitational force
function Fg = gravity(g,m)
    Fg = g*m;
end

% Function fr buoyancy force
function Fb = buoyancy(d,g,l0,rho)
    Fb = g*rho*l0*d^2/4;
end

% Function for total mass of node
function Fm = mass(d,l0,m,rho)
    Fm = m+eye(3)*rho*pi*d^2/4*l0;
end

% Function for the spherical boundaries
function fa = fa_s(s_s,x_s,y_s,z_s,a_s,b_s,c_s)
    fa = ((- a_s^2*y_s^2 - a_s^2*z_s^2 + 2*a_s*b_s*x_s*y_s + 2*a_s*c_s*x_s*z_s - b_s^2*x_s^2 - b_s^2*z_s^2 + 2*b_s*c_s*y_s*z_s - c_s^2*x_s^2 - c_s^2*y_s^2 + s_s^2*x_s^2 + s_s^2*y_s^2 + s_s^2*z_s^2)^(1/2) + a_s*x_s + b_s*y_s + c_s*z_s)/(x_s^2 + y_s^2 + z_s^2);
end
