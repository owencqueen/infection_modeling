# Challenge 19.4

### Prompt
"Run the model of Challenge 19.3 1000 times, recording the number of individuals who become infected in each run. (Note that this is equal to the number of recovered individuals when the run is terminated.) Plot this data as a histogram (as in Figure 19.2), and compute the mean proportion of recovered individuals and the variance in this number. Try several different values of ν to see whether the variance changes."

## Methods
This challenge required a repeated run of the model previously written in Challenge 19.3. In order to do this, I wrote a function that varied the values of ν while repeatedly running the model and recording the results of each run in a Pandas dataframe. We know that the final number for this represents the total number of infected patients throughout the course of the model run because the only pathway out of infection is recovery. Therefore, I plotted the frequencies of the recovered individuals on these histograms.

## Results
I displayed the data from these runs in histograms similar to those provided as an example in the writeup. I added the ν value, the mean proportion of infected individuals, and the variance of the resulting infected proportions from the model. This data is displayed below:

![Varying Nu Values](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-4/infected_histograms.png)

## Discussion
We can clearly see that a larger ν value decreases the spread of our infectious disease in the model. However,  we can see that the reduction in mean proportion of infected patients is nonlinear as ν decreases. This finding will be discussed more in the next challenge. 
