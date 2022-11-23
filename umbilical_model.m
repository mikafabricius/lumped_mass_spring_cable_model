function [r,v,a] = umbilical_model(c_length,segments,v_ship,current,waves,Ts)
    
    % Make all parameters global, to be used in local functions
    global Cn Ct d dt E g k l0 lc m n rho vc mc

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
    t = 0:dt:100; % s

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
    end

    % Initial conditions, important to take into account the secant method,
    % since it actually needs unique initial guesses of r_i,r_i-1 and r_i+1.
    % This means that three different initial guesses must be made in order for
    % the algorithm to converge.
    
    r0 = [0;0;0]; % m, initial position of first node
    v0 = [0;0;0]; % m/s, initial velocity of first node
    a0 = [0;0;0]; % m/s^2, initial acceleration ofr first node
    
    r(1,:,1) = r0;
    v(1,:,1) = v0;
    a(1,:,1) = a0;
    
    % Set the velocity of the first node at every timestep, equal
    % to the ships velocity
    v(1,:,:) = v_ship;
    
    % Initial guess assuming no acceleration or velocity. 
    % Just assuming freefall.
    
    for i=1:length(1:n)
        r(i+1,:,1) = r(i,:,1) + [l0,l0,0];
        r(i+1,:,2) = r(i+1,:,1) - [0.00001,0.00001,0];
        r(i+1,:,3) = r(i+1,:,2) - [0.00001,0.00001,0];
        r(i+1,:,4) = r(i+1,:,3) - [0.00001,0.00001,0];
    end

    % Display a message in the command window that tells the user that the
    % initialization is complete.
    
    msg = ('Initialization complete.');
    disp(msg)

    %% Runge-Kutta fourth order algorithm
    
    % Set timestep to be h, to follow the standard Runge-Kutta notation.
    h = dt;
    
    % Start for loop, to loop through every timestep.
    for i = 3:length(t)-1
        
        % Start second for loop, to loop through every node except node 1
        % and n+1.
        
        for j=2:n
            
            % Predict the position of the first node/the ship using forward
            % Euler method.
            
            r(1,:,i+1) = r(1,:,i)+h*v(1,:,i);

            % Forward Euler method for predicting the i+1 timestep position

            r(j,:,i+1) = r(j,:,i) + h*v(j,:,i);
            
            % Use the system of first order ODEs, from the Runge-Kutta method,
            % in order to correct the prediction and calculate the velocity and
            % acceleration.
            
            K1r = h*f1(t(i),r(j,:,i),v(j,:,i));
            K1v = h*f2(t(i),r(j,:,i),v(j,:,i),r(j,:,i+1),r(j,:,i-1),r(j-1,:,i+1),waves);
            K2r = h*f1(t(i)+h/2,r(j,:,i)+K1r/2,v(j,:,i)+K1v/2);
            K2v = h*f2(t(i)+h/2,r(j,:,i)+K1r/2,v(j,:,i)+K1v/2,r(j,:,i+1),r(j,:,i-1),r(j-1,:,i+1),waves);
            K3r = h*f1(t(i)+h/2,r(j,:,i)+K2r/2,v(j,:,i)+K2v/2);
            K3v = h*f2(t(i)+h/2,r(j,:,i)+K2r/2,v(j,:,i)+K2v/2,r(j,:,i+1),r(j,:,i-1),r(j-1,:,i+1),waves);
            K4r = h*f1(t(i)+h,r(j,:,i)+K3r/2,v(j,:,i)+K3v/2);
            K4v = h*f2(t(i)+h,r(j,:,i)+K3r/2,v(j,:,i)+K3v/2,r(j,:,i+1),r(j,:,i-1),r(j-1,:,i+1),waves);
            
            % Updating the r(j,:) and v(j,:) values using the Runge-Kutte
            % method.
            
            r(j,:,i+1) = r(j,:,i+1) + 1/6*(K1r+2*K2r+2*K3r+K4r);
            v(j,:,i+1) = v(j,:,i) + 1/6*(K1v+2*K2v+2*K3v+K4v);
            a(j,:,i+1) = K1v;
            
            % Print message for debugging purposes
            
            %msg = [num2str(j),'th node done at ',num2str(i),'th timestep.  Simulation time: ',num2str(i*dt),'[s]'];
            %disp(msg)
            
            if abs(buoyancy()-gravity()-drag(v(j,:,i),r(j,:,i+1),r(j,:,i),r(j,:,i-1))) > abs(tension(r(j,:,i),r(j,:,i+1),r(j-1,:,i+1))-tension(r(j,:,i-1),r(j,:,i),r(j-1,:,i+1)))
                
                % Do not apply spherical boundary conditions if the forces
                % applied to the segment, stretches the segment.
                
            elseif norm(r(j-1,:,i+1)-r(j,:,i+1)) > l0 || norm(r(j-1,:,i+1)-r(j,:,i+1)) < l0
                
                % Apply the spherical boundary conditions, since the segment is
                % not supposed to stretch.
                
                P = r(j,:,i+1) - r(j-1,:,i+1);
                Q = (l0/norm(P))*P;
                r(j,:,i+1) = Q + r(j-1,:,i+1);
            end
        end

        msg = ['all nodes evaluated at ',num2str(i),'th timestep.  Simulation time: ',num2str(i*dt),'[s]'];
        disp(msg)
    end
    %% Plot for position of nodes
    for i=1:length(t)-1
        plot3(r(1:end-1,1,i),r(1:end-1,2,i),r(1:end-1,3,i),'*-')
        xlim([-20 20])
        ylim([-20 20])
        zlim([-20 20])
        grid on
        if i == 1
            pause(1)
        else
            pause(0.01)
        end
    end
end
%% Functions for equation of motion and functions for Runge-Kutta method

% Function for the first, 1st order ODE for the RK method

function f1 = f1(t,r_i,v_i)
    f1 = v_i;
end

% Function for the second, 1st order ODE for the RK method

function f2 = f2(t,r_i,v_i,r_ip1,r_im1,r_jm1ip1,waves)
    if abs(buoyancy()-gravity()-drag(v_i,r_ip1,r_i,r_im1)) > abs(tension(r_i,r_ip1,r_jm1ip1)-tension(r_im1,r_i,r_jm1ip1))
        f2 = (buoyancy()-gravity()-drag(v_i,r_ip1,r_i,r_im1)+tension(r_i,r_ip1,r_jm1ip1)-tension(r_im1,r_i,r_jm1ip1))/(mass());
    else
        f2 = (buoyancy()-gravity()-drag(v_i,r_ip1,r_i,r_im1)+wave(waves))/(mass());
    end
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

function Ft = tension(r_i,r_ip1,r_jm1ip1)
    global d E l0
    Ft = pi*d^2*E*(norm(r_ip1-r_i)-l0)*(r_ip1-r_i)/(4*l0*norm(r_ip1-r_i));
    
    % If it is compression that occurs, the tension should b 0
    
    if norm(r_ip1-r_jm1ip1) < l0
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

function Wd = wave(waves)
    Wd = waves;
end