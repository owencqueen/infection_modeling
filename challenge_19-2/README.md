# Challenge 19.2

### Prompt
"Modify your model to include mobility and run it for δ = 0.01 until there are no infected patients. Display the results as in Challenge 19.1."

## Methods
This challenge introduces the concept of mobility, which is the ability for patients to swap positions with other patients on our patient graph. This mobility is represented by δ, which is the probability that a given patient may swap with another patient. 

When I was implementing the mobility feature of the model, the matrix representation of the patient graph introduced a unique problem for the mobility process. My simulation works by iterating over the patient graph, and this happens by indexing the patients. When the patients move, this means their indices change, so there is a chance that this individual patient could undergo a state change more than once. For example, if a patient moves from index [2,2] to index [3,3], this patient is essentially iterated over twice. This would give this patient an unequal probability to change states, which is an undesired result of the model. In addition, a patient could move out of the line of iteration, giving it a 0% chance of undergoing a state change.

To circumvent this issue, I decided to implement the mobility independent of the other simulation processes. During the simulation, we first run the mobility process to completion, and then we proceed to the stochastic process of state changing. This means that the patients are stationary during the iteration, which is the result that we want to achieve in our model. This preserves the probability of each patient undergoing a state change independently.

## Results
The plot generated from the challenge was very similar to the plot from Challenge 19.1. Upon running this model for several trials, it appears as if there is more variation in the results when mobility is implemented versus when it is not implemented (as in Challenge 19.1). However, I have only included a plot from one trial.:

![Mobility Plot](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-2/model_with_mobility.png)

### Testing Uniformity of Random Numbers
In order to preserve the probabilities of state changes in this model, we need to be sure that our random numbers are drawn from a uniform distribution. For this project, I decided to verify the uniformity of choosing my random indices in the function `rand_ind1` (used in the mobility procedure). The plot below shows the result of running my method for 50,000 trials. We can see that our frequency of choosing indices is relatively uniform, so the conclusion was made that my method for deterimining random indices was based on uniform probability distribution.

![Uniform Indices](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-2/plot_uniformity_test.png)

## Run the Script
Running this script requires no additional actions. This script will output a plot of the counts of Infected, Recovered, and Susceptible patients to the screen.
