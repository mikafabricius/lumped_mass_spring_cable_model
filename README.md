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
![image](https://user-images.githubusercontent.com/26135452/202176712-ddb77f4f-4406-43bf-b84a-29e7deb860b9.png)

![image](https://user-images.githubusercontent.com/26135452/202176661-44f57ef3-1b87-4969-850a-a05ae6584e09.png)

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
