# Lumped-mass-spring Umbilical model
Matlab code for an LMS umbilical model, with spherical bound nodes for linkage of segments. Using fourth-order Runge-Kutta and secant method for finding solution numerically.

# MATLAB file umbilical_model.m contains

- Detailed comments for explanation
- Detailed explanation of the equation of motion
- Use and navigation of 3 dimensional matrices in MATLAB
- Forward Euler method for for prediction step
- Fourth-order implicit Runge-Kutta method
- Plots for visualization of umbilical dynamics


Example of plots:
![image](https://user-images.githubusercontent.com/26135452/202176712-ddb77f4f-4406-43bf-b84a-29e7deb860b9.png)

# To use the model:

step 1: Give umbilical parameters.
- Umbilical weight, given in weight per km of umbilical cable.
- Umbilical cable length
- Umbilical cable diameter
- Young's modulus if using tension
- Tangential and normal drag coefficients (Can be estimated through the use of Reynolds number)

step 2: Give model parameters.
- Timestep
- Umbilical nodes/segments
- Simulation runtime

step 3: Run the simulation using the function call:

    [r,v,a] = umbilical_model(cable_length,segments,v_ship,current,waves,Ts)


# Specifically for use in simulink simulation

step 1: Go in simulink and create matlab function block
step 2: Copy the code from the file umbilical_model_simulink.m into the matlab function block in simulink
step 3: Setup standard input parameters in simulink (t,clength,v_ship,current,waves,Ts)
step 4: Change the amount of nodes (variable n), inside the matlab function block
step 5: Setup remainder of inputs in the following ways:
- Setup forward euler method (r_ip1 = r_i + Ts*v_i) for input into r_ip1
- Feedback the function outputs into r_i, r_im1 and v_i
Step 6: Add unit delay blocks in order to break algebraic loops
- Add unit delay into r_i and v_i
- Add two unit delays into r_im1
- The following figure shows an example setup

<img width="540" alt="Screenshot 2022-11-16 at 13 54 44" src="https://user-images.githubusercontent.com/26135452/202186396-75c0345a-a86b-470b-b797-5e6b30524ad1.png">

step 7: Set initial conditions for the unit delay blocks using the following code:

            % Create initial conditions
            r_ini_i = zeros(n,3);
            v_ini_i = zeros(n,3);
            r_ini_im1 = zeros(n,3);

            for i=1:n-1
                r_ini_i(i+1,:) = r_ini_im2(i,:) + [l0,0,0];
                r_ini_im1(i+1,:) = r_ini_im2(i+1,:) - 0.1;
            end

Where n is the amount of nodes and l0 is the initial length of each segment.
This just assumes freefall for the initial conditions and should be adjusted according to your setup.
