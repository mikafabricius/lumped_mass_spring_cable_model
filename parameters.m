clear 
clc 

% Define the parameters of the cable

n = 20; % # of segments
k = int32(n+1); % # of nodes
g = [0,0,9.81]; % m/s^2, gravitational acceleration
rho = 1020; % kg/m^3, density of seawater
d = 0.01; % m, diameter of umbilical
E = 200*10^9; % Pa, kg/(m*s^2), Young's modulus
lc = 20; % m, umbilical cable length
m_km = 250; % kg/km, weight per km umbilical cable
mc = m_km * lc/1000; % kg, total weight for umbilical cable
Ct = 1; % Tangential drag coefficient
Cn = 1.2; % Normal drag coefficient
l0 = lc/n; % initial length of all segments
m = eye(3)*mc/n; % kg, weight per segment

% Create initial conditions
r_ini_i = zeros(k,3);
v_ini_i = zeros(k,3);


for i=1:k-1
    r_ini_i(i+1,:) = r_ini_i(i,:) + [l0,0,0];
end