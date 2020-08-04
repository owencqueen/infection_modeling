# Challenge 19.1

### Prompt:
"Run this model for m = n = 10, k = 4, τ = 0.2 until there are no infected patients. Plot I(t), S(t), and R(t) in a single graph.
If possible, display the epidemic as a movie. To do this, form a matrix of dimension m × n, where the value of each entry corresponds to the state of the corresponding patient on a particular day. Using the movie command, we can display these matrices in sequence, day after day."

## Methods
Since this was the first challenge, I needed to construct my model before I could obtain any results. I decided to take a very object-oriented approach to the model construction, utilizing a class in Python to store all of the parameters for the model as well as functions that would run our simulations. 

### Patient Graph

![Graph of People](https://github.com/owencqueen/infection_modeling/blob/master/imgs/grid_of_people.jpg)

In order to represent my patient, I chose to use a Numpy array instead of constructing a traditional graph. Since our graph contains obvious symmetry and patterns (all adjacent patients are connected to each other), then a graph representation would have introduced unnecessary complexity and storage. 

### Monte Carlo Simulation
This class contains functions that run a Monte Carlo simulation on the patient graph. This simulation proceeds by drawing random numbers and using our infection probability (τ) to decide whether a susceptible patient adjacent to an infected patient becomes infected. This simulation will later be expanded upon in the other challenges.

## Results
In this portion of the project, I only generated one plot as the prompt directed. This plot is shown below, and in order to see the movie (sequence of patient graphs), one can run the script **[run_challenge_19-1.py](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-1/run_challenge_19-1.py)**.

![Challenge 19.1 Plot](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-1/infection_simulation_19-1.png)

## Running the Code
Running this script requires no additional actions. This script will output a plot of the counts of Infected, Recovered, and Susceptible patients to the screen (with mobility included).
