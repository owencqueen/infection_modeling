# Challenge 19.5

### Prompt
"Develop a vaccination strategy that, on average, limits the epidemic to 20% of the population. Do this by using a nonlinear equation solver to solve the problem `R(ν)−.2 = 0`, where R(ν) is the mean proportion of recovered individuals when we use a vaccination rate of ν. For each value of ν presented by the solver, you need to get a reliable estimate of R by running the model multiple times. Use the variance estimates from Challenge 19.4 to determine how many runs to use, and justify your choice."

## Methods/Rationale
This challenge required me to write a method that would find a function to describe our data points for varying ν values. After finding this function, I then needed to find the point at which this function was equal to the percent at which we wanted to contain the outbreak. In order to find this value, the problem statement suggests exploiting root-finding methods that I have learned in COSC 377. 

### Choice of Curve Fitting
One major choice that I had to make for this challenge was the method that I would use to compute a function to describe my data. I considered interpolation, but I wanted to use a streamlined method that would be easy to work with when I was finding roots of the equation. Therefore, I considered techniques that fit polynomial curves. I have recently taken BAS 320 where I learned a lot about regression modeling and least-squares fitting. Therefore, I decided to use least-squared fitting, and I restricted my polynomial to a maximum of degree 9 because any fit over degree 9 might cause extreme extrapolation. 

### Choice of Optimization Parameter
Once I had chosen to use least-squares fitting, I needed to find the optimal degree for the polynomial fit. I needed a polynomial degree that was high enough to fit our nonlinear trend, but I wanted to avoid fitting a function that caused extrapolation outside the bounds of our x values. Therefore, I needed a metric that would determine how well our function fit the data while also penalizing the function for having a high degree. Therefore, I considered three different goodness-of-fit metrics for my model selection: Root mean squared error (RMSE), R^2 (pearson's correlation coefficient squared) adjusted (R2adj), and Akaike information criterion (AIC). 

#### RMSE
The RMSE is a basic measure of goodness-of-fit. This metric simply takes the square root of the sum of squared errors or the difference between the predicted and actual values, also know as the residual. This metric does not penalize the model based on degree of the polynomial, so I used this metric as a control when comparing the results to other metrics that took degree of polynomial into consideration. The model with the lowest RMSE value was consistently the **degree 9** polynomial.

#### R2adj
The R2adj metric is a bit more sophisticated than the RMSE. Based off of the classic goodness-of-fit measure, R^2, the R2adj statistic adjusts the goodness-of-fit test to account for the number of parameters in our model, which in our case is the degree of the polynomial. Because of this factor, the R2adj seemed to be a better measure of goodness-of-fit than RMSE. When using this statistic for comparison, we also consistently chose the **degree 9** polynomial. 

#### AIC
Even more sophisticated than R2Adj is the AIC statistic. This metric is commonly used to choose multiple regression models as well as machine learning models in some instances. AIC also accounts for number of parameters in our model. This metric seemed to be the most accurate in terms of penalizing the model for having a high degree. However, the model with the highest (best) AIC was also consistently the **degree 9** polynomial. 

#### Conclusion
Every goodness-of-fit measure that I used to determine the optimal degree of my polynomial fit showed that the degree 9 polynomial was the best overall model for predicting the data. At first, I thought that this finding invalidated all of the work I had done to implement these different measures; however, there is significant evidence to believe that the degree 9 polynomial is the best choice for calculating our root. 

First, our data has very little noise; i.e. the repeated values tend to correlate strongly in a nonliner relationship. Therefore, our polynomial, while being of a high degree, does not oscillate between values. This phenomenon can be seen in the plot below which shows the degree 9 polynomial fit to 100 evenly spaced nu values (fit using R2adj optimization).:

![Least Squares Curve Fit](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-5/plots/R2adj_fit_optimization_curve.png)

Second, our data presents a unique situation that guards against extrapolation. Our range of nu values is restricted to [0,1], and in this analysis, 100 values are considered in that range, with a minimum value of 0 and maximum value of 1. Therefore, we would not be using this curve to predict values outside of the range on which the curve was fit. This means that the extrapolation would have to occur within this range, but since we already fit our curve within this range using many data points, this makes the chances of significant extrapolation very low. 

Therefore, considering both the two considerations above as well as the validation of three proven goodness-of-fit metrics, I chose to go forward with this challenge using the degree 9 polynomial. 

### Choice of Root-finding techinique
To find the root of this function, I decided to use the Newton-Raphson algorithm. This is a quick, efficient, and accurate algorithm, but normally, its caveat is that it requires one to compute the derivative for the function. However, since I fitted my curve with the least squares technique, my function was a polynomial, so it was easy and efficient to compute its analytic derivative. Therefore, I wrote a function to compute the derivative of my fitted function, and the Newton-Raphson algorithm was easy to implement by using the provided module. 

## Results

### Optimal nu for 20% Containment
After running this model for 20 trials and taking the average optimized nu value (found by the root-finding method), I have found that with our given parameters:

m = 10,
n = 10,
k = 4,
tau = 0.2,
delta = 0.01,
nu = 0.1

A value of nu = 0.0287 would contain the outbreak to 20% of the population of patients.

### Examining nu Values for other Containment Levels
After the above finding, I decided to examine the difference in nu values needed if we were to contain the outbreak to differnt percentages. My findings are shown in the plot below:

![Varying nu for containment](https://github.com/owencqueen/infection_modeling/blob/master/challenge_19-5/plots/nu_for_different_cont_pct.png)

## Discussion
We can clearly see that the optimal nu value decreases in a nonlinear fashion as the percent containment increases. An increase in nu value at low containment percentages, such as in the range [0.03,0.075], results in a small decrease in percent containment for the outbreak. This trend was predicted in the results from Challenge 19.4, but the plot above shows concrete proof of this phenomenon. 

In terms of practical meaning, this finding is significant for public health officials. As we are seeing during the race for a vaccine for COVID-19, mass-producing a vaccine is a very expensive task. Therefore, when we are producing a vaccine, we want to maximize our resources, allocating enough resources for vaccine production but also for production for other medical equipment such as ventillators to treat infected patients. Therefore, we want to produce enough vaccines so that we are effectively defending the population against the infectious disease; however, we also need to make sure that our resources are allocated so that we are producing the minimum amount of vaccines that still makes our defense strategy effective. 

### Herd Immunity
The metric that we have calculated in this challenge could be considered to be the vaccine percentage needed to achieve some level of [herd immunity](https://www.jhsph.edu/covid-19/articles/achieving-herd-immunity-with-covid19.html). Herd immunity is an epidemiological phenomenon that is observed when a sufficient percentage of a population is vaccinated or immune to an infectious disease. What happens at high levels of vaccination is that vaccinated individuals effectively block non-vaccinated individuals from contracting the disease. This phenomenon can be observed in the above graphic, where a vaccination rate of 0.075 results in only 10% of the population contracting the disease. Therefore, a model such as this one could be used to calculate these herd immunity statistics, as we have effectively done in finding our optimal nu for certain percentages of containment.

## Running the Code
Running this script requires some additional actions on the part of the user.

1. Un-comment the include statements in `infection_model.py`
  - In addition to moving `infection_model.py` into the `challenge_19-5` directory, the user needs to remove the comments over lines 12 and 13, which are (as shown in the unedited `infection_model.py` file):
  ```
  # Uncomment these import statements if running Challenge 19.5
  #from polyFit import * 
  #from newtonRaphson import *
  ```
  - This is required to run any of the scripts for Challenge 19.5
  - Once the user is done running this challenge, they should comment these lines back
  
2. Un-comment the desired block of code in the `run_challenge_19-5.py` file
  - There are two main blocks of code in the script file:
   
   **Block #1**: Runs the model and finds the optimal nu value at a percent containment of 0.2
   
   **Block #2**: Runs the model and plots the differing nu values needed for different containment levels
   
  - The user must choose which block to un-comment by removing the block quote (`''' '''` in Python) around the block of code
  
Only after completing these two steps can the user then run this script.
