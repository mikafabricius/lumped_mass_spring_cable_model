# LumpedMassSpring_Umbilical_model
Matlab code for an LMS umbilical model, with spherical bound nodes for linkage of segments. Using fourth-order Runge-Kutta and secant method for finding solution numerically.

# MATLAB file contains

- Detailed comments for explanation
- Detailed explanation of the equation of motion for each node
- Use and navigation of 3 dimensional matrices in MATLAB
- Secant method for prediction step
- Fourth-order Runge-Kutta method
- Plots for visualization of umbilical dynamics


Example of plots:

<img width="522" alt="Screenshot 2022-11-16 at 13 01 15" src="https://user-images.githubusercontent.com/26135452/202175830-054c4815-f80e-49a4-bb25-b47be3dfda3c.png">

# To use the model:

step 1: Give umbilical parameters.
- Umbilical weight, given in weight per km of umbilical cable.
- Umbilical cable length
- Umbilical cable diameter
- Young's modulus (if using tension for boundary conditions)
- Tangential and normal drag coefficients (Can be estimated through the use of Reynolds number)

step 2: Give model parameters.
- Timestep
- Umbilical nodes/segments
- Simulation runtime

% work in progress
