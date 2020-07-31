# Challenge 19.3

### Prompt
"Suppose that each day each susceptible individual has a probability ν of being vaccinated. Rerun your model with ν = 0.1 until there are no infected patients. Display the results as in Challenge 19.1 and compare the results of the three models."

## Methods
This challenge introduced the idea of vaccination, and I decided to implement this within the main Monte Carlo simulation. The vaccine was administered at a probability of ν only to Susceptible patients. 

I experienced a similar problem with the implementation of the vaccination prcoess as in the mobility process: the interation process seemed to introduce unwanted results. If the vaccination process is implemented within the Monte Carlo simulation along with the other state changes, the patients with higher index numbers have a smaller probability of being vaccinated. 

For example, if we are iterating through our patient graph, and we come to patient [5,5] who is infected, patient [5,6] (patient [5,5]'s neighbor) has a probability of being infected if patient [5,6] is susceptible. Therefore, patient [5,6] has a probability of `(τ)*(ν)` of being vaccinated rather than a probability of `(ν)` if the patient proceeded their infected neighbor, such as patient [5,4] in our example (assuming that patient does not have another infected neighbor). This probability is lower which means that the order of iteration matters in the vaccination process. Therefore, I had to consider when I needed to run the vaccination procedure for the patients in our patient graph. I explored this consideration more in Challenge 19.4. 

## Results
The resulting plot from this challenge is shown below:

![Vaccine Plot](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-3/infection_model_vaccine.png)

For convenience, I have included the plots from Challenge 19.1 and 19.2 below, respectively.:

![19.1 plot](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-1/infection_simulation_19-1.png)
#### Challenge 19.1
- This is the basic version of the model
- Excludes both mobility of patients and vaccination of patients


![19.2 plot](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-2/model_with_mobility.png)
#### Challenge 19.2
- Introduces mobility but omits vaccination

## Discussion
Before comparing and contrasting these plots, it is important to note that all of these procedures are stochastic which means that every trial produces a different result. Average results of the model over time are compared in later challenges. 

One metric we could use to measure the extent of the outbreak is the final number of Recovered patients when the outbreak process has finished running. By examining this metric, we can say with reasonable confidence that the model in 19.3 best contains the outbreak, and the model in 19.2 resulted in the worst outbreak. The model in 19.1 resulted in an outcome that was very close to 19.2. 

These results are what we would intuitively expect from the model. The introduction of a vaccination and holding the patients stationary seemed to be effective in limiting the outbreak. These are common techniques used by public health officials during infectious disease outbreaks. On the surface, these models seem to validate these stratgies.

