function [r,v,a] = cable_model_matlab(cableLength,segments,velocityShip,current,waves,ts,simulationTime,r_ini,v_ini)
    %% Define the parameter values
    
    % Values that you should change depending on cable
    
    m_km = 320; % kg/km, weight per km umbilical cable
    E = 200*10^9; % Pa, kg/(m*s^2), Young's modulus
    d = 0.01; % m, diameter of umbilical


    dt = ts; % s, timestep
    vc = current; % m/s, water current
    lc = cableLength; % m, umbilical cable length
    n = segments; % # of segments
    k = n+1; % # of nodes
    g = [0,0,9.81]; % m/s^2, gravitational acceleration
    rho = 1020; % kg/m^3, density of seawater
    mc = m_km * lc/1000; % kg, total weight for umbilical cable
    Ct = 1.2; % Tangential drag coefficient
    Cn = 1.2; % Normal drag coefficient
    l0 = lc/n; % initial length of all segments
    m = eye(3)*mc/n; % kg, weight per segment

    % (Mass + added mass) matrix
    M = m+eye(3)*rho*pi*d^2/4*l0;

    % Simulation time
    t = 0:dt:simulationTime; % s

    %% Generate multidimensional array for storing solution of r,v,a in x,y,z.
    % The third dimension is added as pages, one page for each timestep
    %         x   y   z         x   y   z         x   y   z
    %        ___________       ___________       ___________
    %   r_1 |___|___|___| v_1 |___|___|___| a_1 |___|___|___|
    %   r_2 |___|___|___| v_3 |___|___|___| a_2 |___|___|___|
    %    :  |___|___|___|  :  |___|___|___|  :  |___|___|___|
    %   r_n |___|___|___| v_n |___|___|___| a_n |___|___|___|

    % Initializing three n by 3 matrices
    r = zeros(k,3);
    v = zeros(k,3);
    a = zeros(k,3);

    % Expanding it to a multidimensional array
    for i = 1:length(t)
        r(:,:,i) = zeros(k,3);
        v(:,:,i) = zeros(k,3);
        a(:,:,i) = zeros(k,3);
    end

    %% Initial position and velocity
    r(:,:,1) = r_ini;
    v(:,:,1) = v_ini;

    %% Display a message in the command window that tells the user that the initialization is complete.
    
    msg = ('Initialization complete.');
    disp(msg)
    pause(1)
    %% Fourth order Runge-Kutta fourth order algorithm
    
    % Start for loop to loop through every timestep.
    for i = 1:length(t)-1

        % Use the system of first order ODEs, from the Runge-Kutta method,
        % to find the position, velocity and acceleration at the next
        % timestep

        K1r = dt*f1(v(:,:,i));
        K1v = dt*f2(r(:,:,i),v(:,:,i),waves,vc,Cn,Ct,d,rho,E,l0,g,m,M);
        K2r = dt*f1(v(:,:,i)+K1v/2);
        K2v = dt*f2(r(:,:,i)+K1r/2,v(:,:,i)+K1v/2,waves,vc,Cn,Ct,d,rho,E,l0,g,m,M);
        K3r = dt*f1(v(:,:,i)+K2v/2);
        K3v = dt*f2(r(:,:,i)+K2r/2,v(:,:,i)+K2v/2,waves,vc,Cn,Ct,d,rho,E,l0,g,m,M);
        K4r = dt*f1(v(:,:,i)+K3v/2);
        K4v = dt*f2(r(:,:,i)+K3r/2,v(:,:,i)+K3v/2,waves,vc,Cn,Ct,d,rho,E,l0,g,m,M);
        
        r(:,:,i+1) = r(:,:,i) + 1/6*(K1r+2*K2r+2*K3r+K4r);
        v(:,:,i+1) = v(:,:,i) + 1/6*(K1v+2*K2v+2*K3v+K4v);
        a(:,:,i+1) = K1v;

        % For the first node, use the ship velocity to calculate next
        % position, using forward Euler method for simplicity.
        r(1,:,i+1) = r(1,:,i)+dt*velocityShip;

        %% Insert the last node dynamics here

        %% Display simulation runtime message   
        msg = ['all nodes evaluated at ',num2str(i),'th timestep.  Simulation time: ',num2str(i*dt),'[s]'];
        disp(msg)
    end

    % Display message that the simulation is done
    msg = 'Simulation complete';
    disp(msg)

end


%% Functions for equation of motion and for Runge-Kutta method
% Functions include:
%   Gravity
%   Buoyancy
%   Tension
%   Drag

% Function for the first, 1st order ODE for the 4RK method
function f1 = f1(v)
    f1 = v;
end

% Function for the second, 1st order ODE for the 4RK method
function f2 = f2(r,v,waves,vc,Cn,Ct,d,rho,E,l0,g,m,M)
    f2 = (buoyancy(r,d,g,l0,rho)-gravity(r,g,m)-drag(r,v,vc,Cn,Ct,d,rho)+tension_im1(r,d,E,l0)+tension_i(r,d,E,l0)-waves)/(M);
end

% Function for drag force, using the tangential and normal components
function Fd = drag(r,dr,vc,Cn,Ct,d,rho)
   Fd = zeros(length(r),3);
   for i=1:length(r)
        if i == 1
            Fd(i,:) = [0,0,0];
        else
            if i == length(r)
                t = (r(i,:)-r(i-1,:))/norm(r(i,:)-r(i-1,:));
            else
                t = (r(i+1,:)-r(i,:))/norm(r(i+1,:)-r(i,:));
            end
            D_t = 1/2*rho*Ct*norm(r(i-1,:)+r(i,:))*d*((dr(i,:)-vc).*t).^2;
            D_n = pi/8*rho*d^2*Cn*((dr(i,:)-vc)-(dr(i,:)-vc).*t).^2;
            Fd(i,:) = D_t + D_n;
            Fd(i,:) = Fd(i,:) .* sign(dr(i,:));
        end
    end
end 

% Function for tension force in i-1'th segment
function Ft = tension_im1(r,d,E,l0)
    Ft = zeros(length(r),3);
    for i=1:length(r)
        if i == 1
            Ft(i,:) = [0,0,0];
        else
            Ft(i,:) = pi*d^2/4*E*(norm(r(i-1,:)-r(i,:))-l0)/l0*(r(i-1,:)-r(i,:))/(norm(r(i-1,:)-r(i,:)));
        end
    end
end

% Function for tension force in i'th segment
function Ft = tension_i(r,d,E,l0)
    Ft = zeros(length(r),3);
    for i=1:length(r)
        if i == 1 || i == length(r)
            Ft(i,:) = [0,0,0];
        else
            Ft(i,:) = pi*d^2/4*E*(norm(r(i+1,:)-r(i,:))-l0)/l0*(r(i+1,:)-r(i,:))/(norm(r(i+1,:)-r(i,:)));
        end
    end
end

% Function for gravitational force
function Fg = gravity(r,g,m)
    Fg = zeros(length(r),3);
    for i=1:length(r)
        Fg(i,:) = g*m;
    end
end

% Function fr buoyancy force
function Fb = buoyancy(r,d,g,l0,rho)
    Fb = zeros(length(r),3);
    for i=1:length(r)
        Fb(i,:) = g*rho*l0*d^2/4;
    end
end
