function [r,v,a] = umbilical_model(c_length,segments,v_ship,current,waves,Ts)
    
    % Make all parameters global, to be used in local functions
    global Cn Ct d dt E g k l0 lc m n rho vc waves

    % Load standard variables if none given at function call
    if ~exist(['length','segments','current','waves','Ts'],'var')
        c_length = 20; % m, standard umbilical cable length
        segments = 5; % standard number of segments
        v_ship = [0 0 0]; % m/s, velocity of the ship
        current = [0 0 0]; % m/s, standard current velocities
        waves = [0 0 0]; % given as magnitude of sinusoidal wave
        Ts = 0.1; % s, timestep
    end


    % Define the parameter values
    dt = Ts; % s, timestep
    n = segments; % # of segments
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

    % Simulation time
    t = 0:dt:30; % s

    % Generate multidimensional array for storing solution of r,v,a in x,y,z.
    % The 3 dimension is added as pages, one for each timestep
    %         x   y   z         x   y   z         x   y   z
    %        ___________       ___________       ___________
    %   r_1 |___|___|___| v_1 |___|___|___| a_1 |___|___|___|
    %   r_2 |___|___|___| v_3 |___|___|___| a_2 |___|___|___|
    %   r_n |___|___|___| v_n |___|___|___| a_n |___|___|___|

    % Initializing three 3,3 matrix
    r = zeros(n,3);
    v = zeros(n,3);
    a = zeros(n,3);

    % Expanding it to a multidimensional array
    for i = 2:length(t)-1
        r(:,:,i) = zeros(n,3);
        v(:,:,i) = zeros(n,3);
        a(:,:,i) = zeros(n,3);
        
        % Making the first nodes velocity equal to the ships velocity
        v(1,:,i) = v_ship;
    end

    % Initial conditions, important to take into account the secant method,
    % since it actually needs unique initial guesses of r_i,r_i-1 and r_i+1.
    % This means that three different initial guesses must be made in order for
    % the algorithm to converge.
    r0 = [0;0;0]; % m, position of first node
    v0 = [0;0;0]; % m/s, velocity
    a0 = [0;0;0]; % m/s^2, acceleration

    % Initial guess assuming no acceleration or velocity, and assuming that no
    % nodes have fixed diameter. Just assuming that everything falls in value
    % by 0.1.
    for i=1:length(1:n)
        r(i+1,:,1) = r(i,:,1) + [l0,0,0];
        r(i+1,:,2) = r(i+1,:,1) - 0.1;
        r(i+1,:,3) = r(i+1,:,2) - 0.1;
        r(i+1,:,4) = r(i+1,:,3) - 0.1;
    end
    v(1,:,1) = v0;
    a(1,:,1) = a0;

    % Secant method tolerance
    tol = 1e-3;

    % Spherical boundaries for nodes
    syms s_s k_s x_s y_s z_s a_s b_s c_s
    eq = s_s^2 == (k_s*x_s-a_s)^2 + (k_s*y_s-b_s)^2 + (k_s*z_s-c_s)^2;
    fa_S = solve(eq,k_s);
    fa_S = fa_S(2);

    % Terms in equation of motion
    % G = @() g*m; % Equation for gravity
    % B = @() g*rho*l0*d^2/4; % Equation for buoyancy
    % D = @(v_i,r_ip1,r_i) 2*rho*(v_i-vc).^2*d*(Ct*pi/4+Cn*norm(r_ip1-r_i)); % Equation for drag   
    % T = @(r_ip1,r_i,as) comp*pi*d^2*E*(norm(r_ip1-r_i)-l0)/(4*l0); % Equation for tension
    % M = @() m+eye(3)*rho*pi*d^2/4*l0; % Equation for total mass of segment
    % W_d = @(t) [2,0,0].*[sin(t+pi),sin(t+pi*1/3),sin(t+pi*2/3)]; % Wave disturbance

    msg = ('Initialization complete.');
    disp(msg)

    %% Runge-Kutta fourth order algorithm
    h = dt;
    for i = 3:length(t)-1
        for j=2:n
            % Secant method 
            dx = inf;
            iter = 0;

            % Secant method for predicting the r(j,:,i+1) value, with an
            % initial guess of x2 and x1.
            x1 = r(j,:,i-1);
            x2 = r(j,:,i);% + [0,0,-1];
            % Calculate the r_(i+1) function in the Runge-Kutta method at
            % timestep i-1, this means evaluating all the Delta functions.
            K1r = h*f1(t(i-1),x1,v(j,:,i-1));
            K1v = h*f2(t(i-1),x1,v(j,:,i-1),r(j,:,i),r(j,:,i-2));
            K2r = h*f1(t(i-1)+h/2,x1+K1r/2,v(j,:,i-1)+K1v/2);
            K2v = h*f2(t(i-1)+h/2,x1+K1r/2,v(j,:,i-1)+K1v/2,r(j,:,i),r(j,:,i-2));
            K3r = h*f1(t(i-1)+h/2,x1+K2r/2,v(j,:,i-1)+K2v/2);
            K3v = h*f2(t(i-1)+h/2,x1+K2r/2,v(j,:,i-1)+K2v/2,r(j,:,i),r(j,:,i-2));
            K4r = h*f1(t(i-1)+h,x1+K3r/2,v(j,:,i-1)+K3v/2);
            K4v = h*f2(t(i-1)+h,x1+K3r/2,v(j,:,i-1)+K3v/2,r(j,:,i),r(j,:,i-2));
            %F1 = r_ip1(f1,f2,h,t(i-1),x1,v(j,:,i-1),x1,r(j,:,i-2));
            F1 = 1/6*(K1r+2*K2r+2*K3r+K4r);
            %F1 = 1/6*(f1(0,0,v(j,:,i)) + 2 * (h*f1(0,0,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),r(j,:,i+1),r(j,:,i-1)))/2)) + 2*(h*f1(0,0,v(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(0,0,v(j,:,i)))/2,v(j,:,i)+()/2,r(j,:,i+1),r(j,:,i-1)))/2)))
            while abs(dx) < tol

                % Count iterations
                iter = iter + 1;

                % Evaluate r_(i+1) function in the Runge-Kutta method at x2
                F2 = r_ip1(f1,f2,h,t(i),r(j,:,i),v(j,:,i),x2,x1);
                %1/6*(h*f1(t(i),r(j,:,i),v(j,:,i))+2*h*f1(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2)+2*(h*f1(t(i)+h/2,r(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2,x2))/2,v(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2,x2))/2))+h*f1(t(i)+h,r(j,:,i)+(h*f1(t(i)+h/2,r(j,:,i)+(h*f1(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2))/2,v(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2,x2))/2))/2,v(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2))/2,v(j,:,i)+(h*f2(t(i)+h/2,r(j,:,i)+(h*f1(t(i),r(j,:,i),v(j,:,i)))/2,v(j,:,i)+(h*f2(t(i),r(j,:,i),v(j,:,i),x2))/2,x2))/2,x2))/2));

                % Calculate change in r using secant method
                dx = (F2 - F1)/ (x2 - x1);

                % Second point is now the new first point
                x1 = x2;
                F1 = F2;

                % Update position of second point using half the change in r
                x2 = x2 - dx/2;
            end

            r(j,:,i+1) = x2;
            K1r = h*f1(t(i),r(j,:,i),v(j,:,i));
            K1v = h*f2(t(i),r(j,:,i),v(j,:,i),r(j,:,i+1),r(j,:,i-1));
            K2r = h*f1(t(i)+h/2,r(j,:,i)+K1r/2,v(j,:,i)+K1v/2);
            K2v = h*f2(t(i)+h/2,r(j,:,i)+K1r/2,v(j,:,i)+K1v/2,r(j,:,i+1),r(j,:,i-1));
            K3r = h*f1(t(i)+h/2,r(j,:,i)+K2r/2,v(j,:,i)+K2v/2);
            K3v = h*f2(t(i)+h/2,r(j,:,i)+K2r/2,v(j,:,i)+K2v/2,r(j,:,i+1),r(j,:,i-1));
            K4r = h*f1(t(i)+h,r(j,:,i)+K3r/2,v(j,:,i)+K3v/2);
            K4v = h*f2(t(i)+h,r(j,:,i)+K3r/2,v(j,:,i)+K3v/2,r(j,:,i+1),r(j,:,i-1));
            % Updating the r(j,:,i+1) value using the Runge-Kutta method
            r(j,:,i+1) = r(j,:,i) + 1/6*(K1r+2*K2r+2*K3r+K4r);
            v(j,:,i+1) = v(j,:,i) + 1/6*(K1v+2*K2v+2*K3v+K4v);
            a(j,:,i+1) = K1v;
            % Print message for debugging purposes
            %msg = [num2str(j),'th node done at ',num2str(i),'th timestep.  Simulation time: ',num2str(i*dt),'[s]'];
            %disp(msg)
            if norm(r(j-1,:,i+1)-r(j,:,i+1)) < l0 || norm(r(j-1,:,i+1)-r(j,:,i+1)) > l0
                fa = subs(fa_S,{s_s,x_s, y_s, z_s, a_s, b_s, c_s},{l0,r(j,1,i+1),r(j,2,i+1),r(j,3,i+1),r(j-1,1,i+1),r(j-1,2,i+1),r(j-1,3,i+1)});%k_S(l0,r(j,1,i+1),r(j,2,i+1),r(j,3,i+1),r(j-1,2,i+1),r(j-1,3,i+1),r(j-1,3,i+1))
                %fa = l0 / sqrt((r(j,1,i+1)-r(j-1,1,i+1))^2+(r(j,2,i+1)-r(j-1,2,i+1))^2+(r(j,3,i+1)-r(j-1,3,i+1))^2);
                r(j,:,i+1) = r(j,:,i+1)*fa(1);
            end
        end

        msg = ['all nodes done at ',num2str(i),'th timestep.  Simulation time: ',num2str(i*dt),'[s]'];
        disp(msg)
    end
    %% Plot for position of nodes
    for i=1:length(t)-1
        plot3(r(1:end-1,1,i),r(1:end-1,2,i),r(1:end-1,3,i),'*-')
        xlim([-20 20])
        ylim([-20 20])
        zlim([-20 20])
        grid on
        pause(0.01)
    end
end
%% Functions for equation of motion and functions for Runge-Kutta method

% Function for the first, 1st order ODE for the RK method
function f1 = f1(t,r_i,v_i)
    f1 = v_i;
end

% Function for the second, 1st order ODE for the RK method
function f2 = f2(t,r_i,v_i,r_ip1,r_im1)
    f2 = (buoyancy()-gravity()-drag(v_i,r_ip1,r_i,r_im1)+wave(t))/(mass());
end

% Function for drag force, using the tangential and normal components
function Fd = drag(v_i,r_ip1,r_i,r_im1)
    global Cn Ct d rho vc
    t_i = (r_ip1-r_im1)/norm(r_ip1-r_im1);
    D_t = rho*pi*d/4*Ct*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc).*t_i).*(v_i-vc).*t_i;
    D_n = rho*pi*d/4*Cn*norm((r_i-r_im1)+(r_ip1-r_i))*norm((v_i-vc)-(v_i-vc).*t_i)*((v_i-vc)-(v_i-vc).*t_i);
    Fd = D_t + D_n;
end 

% Function for tension force, currently not used
function Ft = tension(r_i,r_ip1)
    global d E l0
    Ft = pi*d^2*E*(norm(r_ip1-r_i)-l0)*(r_ip1-r_i)/(4*l0*norm(r_ip1-r_i));
    if norm(r_ip1-r_i) < l0
        Ft = 0; 
    end
end

% Function for gravitational force
function Fg = gravity()
    global g m
    Fg = g*m;
end

% Function fr buoyancy force
function Fb = buoyancy()
    global d g l0 rho
    Fb = g*rho*l0*d^2/4;
end

% Function for total mass of node
function Fm = mass()
    global d l0 m rho
    Fm = m+eye(3)*rho*pi*d^2/4*l0;
end

% Function for wave disturbances
function Wd = wave(t)
    global waves
    Wd = waves.*[sin(t+pi),sin(t+pi*1/3),sin(t+pi*2/3)];
end