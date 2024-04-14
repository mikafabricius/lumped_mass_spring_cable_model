# Latest update (14/04/2024)

- (21/12/2022) Fixed a sign error in the function for calculating the drag force. This force led to accumulation of tension forces.
- (14/04/2024) Uploaded a MATLAB GUI for generating initial conditions.


# Lumped-mass-spring cable model

The model is a first principle model, using Newton's second law as the equation of motion:


<img width="564" alt="Screenshot 2022-12-19 at 12 10 51" src="https://user-images.githubusercontent.com/26135452/208413180-4417e290-a58a-473a-b504-f7bcd00af70b.png">


# Simulink simulation file

The Simulink simulation file uses the continuous integrator block in order to integrate the equations of motion, from Newton's second law.
This results in the velocity of each node and the position of each node.


Procedure for getting simulation to work:

step 1: Initialize all variables
- This is quickly done by opening and changing them in the parameters.m file

step 2: Select the correct solver

- Since the model uses tension as the constraint reaction forces, the equation of motion is a stiff equation.
The model has been tested so far with the following solvers, where the fastest solver is the fixed-step solver:
  - Variable-step solvers: ode15s, ode23s, ode23t, ode23tb.
  - Fixed-step solver: ode4 (using ts =< 0.0001)
  
step 3: Set ship velocity, waves and current
  - Ship velocity is set in the parameters.m file
  - Current and waves are set in the simulink file
      
![image](https://user-images.githubusercontent.com/26135452/205927399-02e8ceb7-cd4f-4f6b-a572-f0e928faaf8e.png)


# MATLAB file contains

- Detailed comments for explanation
- Detailed explanation of the equation of motion
- Use and navigation of 3 dimensional matrices in MATLAB
- Fourth-order explicit Runge-Kutta method

# To use the matlab model:

step 1: Set cable parameters
- These parameters are set inside the matlab file at the top of the script

step 2: Insert the last node dynamics
- Insert the equation of motion of the object at the end of the cable at line 86.
  - Depending on what you put at the end of the cable, the last node should have the dynamics of that object.
    If, as an example, the last node is a ROV, insert the equation of motion for the ROV, but remember to use the last tension force as a disturbance in
    the ROV equation of motion.

step 3: Run the simulation using the function call:

    [r,v,a] = cable_model_matlab(cableLength,segments,velocityShip,current,waves,ts,simulationTime,initialPosition,initialVelocity);

- Remember that for the fourth-order Rung-Kutta method the timestep should be ts <= 0.0001

# Plots where the model has been used in an underwater project

- Results of a simulation where the ship, cable and mass are initialized to be horizontally in line with each other and with no velocity from the ship or actuation from the mass. This will cause the negatively buoyant cable to pull the mass down, but since the mass is made positively buoyant, it will stay at the water surface, but instead be pulled inwards, towards the ship.

<img width="1048" alt="Screenshot 2022-12-19 at 12 20 05" src="https://user-images.githubusercontent.com/26135452/208414835-8569f016-0ee9-484e-b5fd-bf2921f9576a.png">

- Results of how the cable interacts with the ship, when the ship is moving at a velocity of 2 m/s in the positive north direction. The cable is 20 meters long, hanging straight down, and is segmented into 10 segments, with 11 nodes. The last node is replaced by a dead mass without any actuation.

<img width="1019" alt="Screenshot 2022-12-19 at 12 20 27" src="https://user-images.githubusercontent.com/26135452/208414872-20807331-7b3a-44a8-ab51-cdfd20bff454.png">
