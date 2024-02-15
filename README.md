# Dynamic Modeler

## Installation
First, clone the repository to a local folder of your choosing.

Package requirements for installation can be found in requirements.txt or environment.yml.

To install requirements using pip, create a new virtual environment. Then run "pip install -r requirements.txt" with the virtual environment active and in the directory of the cloned repository.

## Overview
See [writeup](https://gaoadam.github.io/docs/projects/dynamic_modeler/dynamic_modeler.html) for a more illustrative explanation of the module.

See other [writeup](https://gaoadam.github.io/docs/projects/filters_and_neural_networks/filters_and_neural_networks.html) to see how I use **neural networks** to predict dynamical system behavior using **TensorFlow**.

See generate_signals.ipynb for an example use of the dynamicmodel.py module

The dynamicmodel.py module allows you to simulate a dynamical system and generate signals from it. 

You may want to dig into the dynamicmodel.py module itself to see how dynamical systems are constructed from the x_iterate function. A good example is the x_driven function, which calls the x_iterate function.

## Usage

You may elect to either use your own set of custom "dynamical" functions or use a prebuilt function.

If you elect to use your own system of functions:

* Deftermine your system variables.
* Determine the number of system variables you wish to use. Let's call this number M.
* Define M functions to calculate the system variables' time derivatives as functions of the other variables
* Define an initial state vector of length M, containing the initial values of your variables. There are some conditions for these functions:
    * They must take in their respective system variable
    * They must take in the time value t
    * They must take in a dictionary of custom arguments, called "args". This dictionary can be empty if not needed.
* Pass the following into x_iterate:
    * Initial state vector
    * The value of the timestep
    * The number of timesteps the simulation will iterate through
    * A list containing the M functions for the time derivatives
    * A dictionary containing the custom arguments required for the system variable functions

Again, reference the x_driven function in dynamicmodel.py to see how this is done.

## Images

The following are some plots and "phase portraits" I plotted from my simulation engine.

First I simulate an oscillator with a sinusoidal driving force and a strong damping component:

![damped_oscillator.png](example_plots\damped_oscillator.png)

![damped_oscillator.png](example_plots\damped_phase_portrait.png)

I also simulate an oscillator with multiple sinusoidal driving forces and a negligible damping component:

![damped_oscillator.png](example_plots\damped_oscillator2.png)

![damped_oscillator.png](example_plots\damped_phase_portrait2.png)