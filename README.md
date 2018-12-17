# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---
## OVerview
This project is to build a MPC controller and make the vehicle drive successfully around the track without going out of the track.

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.


## Tips

1. The MPC is recommended to be tested on examples to see if implementation behaves as desired. One possible example
is the vehicle offset of a straight line (reference). If the MPC implementation is correct, it tracks the reference line after some timesteps(not too many).
2. The `lake_track_waypoints.csv` file has waypoints of the lake track. This could fit polynomials and points and see of how well your model tracks curve. NOTE: This file might be not completely in sync with the simulator so your solution should NOT depend on it.
3. For visualization this C++ [matplotlib wrapper](https://github.com/lava/matplotlib-cpp) could be helpful.)
4.  Tips for setting up your environment are available [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)

[//]: # (Image References)

[image1]: ./img/equations.png "Equations"

## Implementation
The MPC controller is implemented in [./src/MPC.cpp](./src/PID.cpp). It utilizes the IPOPT and CppAD libraries to calculate the optimal trajectory and the corresponding actuation commands i.e. the throttle/brake and steering angle, to minimize the cost function of cross track error, steering angle error and penalization of roughness .

### The Model
The MPC model uses a kinematic model without taking into account of the complex road and tire interactions. The model equations are below:
![Equations][image1]
The model has 6 state variables:
- x : car's position x
- y : car's position y
- psi : car's heading angle
- v : car's velocity
- cte : cross track error
- epsi : heading angle error
The variable Lf is the length from front to center of gravity, it is given by Udacity sample code.
The mode has two outputs:
- a : acceleration/deceleration value
- delta: steering angle

The objective is to find the best a and delta values to minimize the cost function of multiple factors:
- cross track error and heading angle error
- penalization of uses of actuations
- penalization of rapid changes

### Timestep Length and Elapsed Duration (N & dt)
Number of points (N) and time interval (dt) together define the prediction horizon. With too many points, the model will become slower easily. However, too less points will not predict a good curve. This project uses the Udacity suggested values: N = 10; dt = 0.1;

### Polynomial Fitting and MPC Preprocessing
In this project, the provided waypoints are transformed to the vehicle coordinates first, then fitted a 3rd order polynomial.

### Model Predictive Control with Latency
To take care of the actuator latency, this project calculates and uses the delayed states to feed into the MPC solver, instead of using the initial values.

## Simulation
Here is a video shows the car successfully drove one lap using the implemented MPC controller: [./videos/mpc.mp4](./videos/mpc.mp4)
