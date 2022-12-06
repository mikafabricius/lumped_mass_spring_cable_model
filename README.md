# Lumped-mass-spring umbilical model

The model is a first principle model, using Newton's second law as the equation of motion:

<img width="479" alt="Screenshot 2022-12-05 at 13 08 56" src="https://user-images.githubusercontent.com/26135452/205633942-48ccddf3-6d07-4dfe-8acf-99eccda88e8a.png">

Currently the MATLAB file is a work in progress...

# Simulink simulation file

The Simullink simulation file uses the continuous integrator block in order to integrate the equations of motion, from Newton's second law.
This results in the velocity of each node and the position of each node.


Procedure for getting simulation to work:

step 1: Initialize all variables
- This is quickly done by opening the "parameters.m" file
step 2: Select the correct solver
- Since the model uses tension as the constraint reaction forces, the equation of motion is a stiff equation.
The model has been tested so far with the following solvers:
Type: Variable-step. Solvers: ode15s, ode23s, ode23t, ode23tb.
Other solvers have not been tested.
      
![image](https://user-images.githubusercontent.com/26135452/205927399-02e8ceb7-cd4f-4f6b-a572-f0e928faaf8e.png)


# MATLAB file contains

- Detailed comments for explanation
- Detailed explanation of the equation of motion
- Use and navigation of 3 dimensional matrices in MATLAB
- Forward Euler method for for prediction step
- Fourth-order implicit Runge-Kutta method
- Plots for visualization of umbilical dynamics

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


