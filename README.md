# Lumped-mass-spring umbilical model

The model is a first principle model, using Newton's second law as the equation of motion:

<img width="479" alt="Screenshot 2022-12-05 at 13 08 56" src="https://user-images.githubusercontent.com/26135452/205633942-48ccddf3-6d07-4dfe-8acf-99eccda88e8a.png">

Currently the MATLAB file is a work in progress...

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


# Simulink simulation file

The Simullink simulation file uses the continuous integrator block in order to integrate the equations of motion, from Newton's second law.
This results in the velocity of each node and the position of each node.


Procedure for getting simulation to work:

step 1: Initialize all variables
    - This is quickly done by opening the "parameters.m" file
step 2: Select the correct solver
    - Since the model uses tension as the constraint reaction forces, the equation of motion is a stiff equation.
      The model has been tested so far with the following solvers:
            Type: Variable-step| Solvers: ode15s, ode23s, ode23t, ode23tb
      The model is numerically stable using some other solvers, like backward-euler for fixed step. The model does however NOT behave as intended, for a       timestep of size 0.1 or larger. <-- Will update on this
      
Figure comming soon.
