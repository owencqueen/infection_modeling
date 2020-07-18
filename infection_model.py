#!/usr/bin/env python3

import math
import numpy as np
import pandas as pd

from collections import Counter
from matplotlib import legend
from matplotlib import pyplot as plt


def rand_ind1(self):
    '''
        model = infection_model()
        x, y = model.rand_ind1()

        Generates random indices within the patient graph of an 
          infection_model object (based on given parameters in model).
        
        Defined outside of class to allow setting of different functions
          to generate random indices.
    '''

    np.random.seed()

    # Find the random swap coordinates
    x = np.random.randint(0, self.n)
    y = np.random.randint(0, self.m)

    return x, y


class infection_model:
    '''
        model = infection_model(m, n, k, tau, delta, nu)

        Generates a model with given parameters

    '''

    random_indices = rand_ind1 # Set random indices choice fn


    def __init__(self, m, n, k, tau, delta = 0, nu = 0):
        # Class constructor - see set_parameters for more information
        self.set_parameters(m, n, k, tau, delta, nu)


    def set_parameters(self, m = False, n = False, k = False, 
                        tau = False, delta = False, nu = False):
        '''
        - Allows user to set/reset parameters after object has been created.
        - Automically resets the internals of the object instance when called.
        - User can give only one parameter at a time. e.g.:
            model = infection_model(m, n, k, tau, delta, nu)
            model.set_parameters(m = 3)
        '''

        if (m != False):
            self.m = m

        if (n != False):
            self.n = n

        if (k != False):
            self.k = k

        if (k != False):
            self.tau = tau
        
        if (delta == False):
            self.mobility = False
        else:
            self.delta = delta
            self.mobility = True

        if (nu == False): # Sets vaccination rate
            self.vaccinate = False
        else:
            self.nu = nu
            self.vaccinate = True

        self.reset_internal()

    # Runs the Monte Carlo simulation for given parameters
    def __monte_carlo_sim(self, prt = False):
        '''
        Private method of class.

        Usage:
        prt: specifies whether to print the patient graph after each 
        iteration.
            - True: print graph
            - False: do not print graph

        Runs a simulation with given parameters over the patient graph.

        - Continues until there are no more infected patients.
        - Updates internal counters of infected, susceptible, recovered,
          and vaccinated patients. 
        '''

        np.random.seed()

        days = 0

        if(prt):
            print(self.patient_graph)

        while True:

            if self.mobility:
                self.__move_patients()

            # Implements infection
            for i in range(0, self.n):
                for j in range(0, self.m):
                    
                    # Case: patient has recovered from infection - R
                    if (self.patient_graph[i][j] == -1):
                        continue

                    # Case: patient if susceptible to infection - S
                    elif (self.patient_graph[i][j] == 0):
                        if (self.vaccinate):

                            vac_prob = np.random.rand()
                            if (vac_prob < self.nu):
                                self.patient_graph[i][j] = -2
                                # Set number to -2 if the patient is vaccinated

                                # Note: has potential for a patient that is later in the
                                # sequence of evaluation to have a lower chance of being 
                                # vaccinated than one earlier in the sequence. Could be 
                                # fixed by (1) allowing infected persons to be vaccinated
                                # or (2) by doing an independent vaccination before this
                                # process.

                    # Case: patient is in disease progression - I
                    elif (self.patient_graph[i][j] > 0):
                        self.patient_graph[i][j] += 1

                        # Generate neighbor indices:
                        neighbors = self.__find_neighbors(i, j)

                        for n in neighbors:
                            # Need to decide whether to infect if neighbor susceptible
                            if (self.patient_graph[ n[0] ][ n[1] ] == 0):
                                rf = np.random.rand() # rf == "random float"
                                # rand() picks rf in range 0-1 by default

                                if (rf < self.tau):   # Stochastic infection of the neighbor
                                    self.patient_graph[ n[0] ][ n[1] ] = 1 

                        # Sub-case: patient done with disease
                        if (self.patient_graph[i][j] == self.k):
                            self.patient_graph[i][j] = -1
                            
            if (prt): # Print graph if specified by user
                print(self.patient_graph)

            self.__update_ISR() # Updates internal counters

            if (self.It[-1] == 0): # If no more infected individuals, break
                return days
            else:
                days += 1
    

    def __update_ISR(self):
        # Updates the I, S, and R counts

        # Flatten the patient graph, count its elements
        patient_count = Counter(self.patient_graph.flatten())

        total_patients = (self.m * self.n)
        
        # Given patient codes: R == -1, S == 0, I == other
        current_R_count = patient_count[-1]
        current_S_count = patient_count[0]

        all_others = current_R_count + current_S_count
            # Counts all the other patients that aren't I

        if self.vaccinate:
            current_V_count = patient_count[-2]
            all_others += current_V_count

        current_I_count = total_patients - all_others
            # We know that if they aren't S or R (or maybe V), they're I

        # Append the counts to our lists
        self.It.append(current_I_count / total_patients)
        self.St.append(current_S_count / total_patients)
        self.Rt.append(current_R_count / total_patients)

        if self.vaccinate:
            self.Vt.append(current_V_count / total_patients)

    
    def __find_neighbors(self, i, j):
        # Finds the indices for the neighbors around one patient
        # i: row index
        # j: col index
        
        neighbors = [(i + 1, j), (i, j + 1), (i - 1, j), (i, j - 1)]

        if (i == 0):
            neighbors.remove((i - 1, j))

        elif ((i + 1) == self.m):
            neighbors.remove((i + 1, j))

        if (j == 0):
            neighbors.remove((i, j - 1))

        elif((j + 1) == self.n):
            neighbors.remove((i, j + 1))

        return neighbors

    def __move_patients(self):
        # Swaps patients at once
        #   Swapped patients chosen stochastically

        np.random.seed()

        moves_dict = {} # Dictionary will hold new locations of each swap

        # For every patient in the patient graph, decide if swap
        # Need to do this independently of actual swapping
        for i in range(0, self.n):
            for j in range(0, self.m):
                rf = np.random.rand()
                if (rf < self.delta): # Moving procedure

                    x, y = self.random_indices()

                    # Initiate the swap - i.e. add to dictionary
                    moves_dict[tuple((i, j))] = tuple((x, y))

        # Swap each patient that needs to be swapped
        while (len(moves_dict) > 0):
            # Choose element in to_move
            #   - Choice of element is arbitrary
            move_from = list(moves_dict.keys())[0]
            
            move_to = moves_dict[ move_from ] # Get place to move to

            # swap( move_from, move_to ):
            temp_copy = self.patient_graph[ move_from[0] ][ move_from[1] ]
            
            self.patient_graph[ move_from[0] ][ move_from[1] ] = \
                self.patient_graph[ move_to[0] ][ move_to[1] ]
            
            self.patient_graph[ move_to[0] ][ move_to[1] ] = temp_copy

            print("moved", move_from, "to", move_to)

            del moves_dict[move_from] # Remove pair we just swapped

            # if tuple_2 in keys: remove that key pair, add new
            if move_to in moves_dict.keys():
                val = moves_dict[move_to] # Get that key-value pair's value:
                
                del moves_dict[move_to] # Remove old k-v pair 

                # Note: move_from is move_to's new location
                moves_dict[move_from] = val # Add new k-v pair

    def __vaccinate_patients(self):
        # Vaccinates patients at once 
        #   - Vaccinated patients chosen stochastically

        # Go through all patients
        # Perform vaccinations independently
        for i in range(0, self.n):
            for j in range(0, self.m):
                vac_prob = np.random.rand()

                if (vac_prob < self.nu):
                    self.patient_graph[i][j] = -2
                    # Vaccination code is -2

    def reset_internal(self):
        # Initialize the patient graph
        self.patient_graph = np.zeros((self.n, self.m))

        # Set patient-zero
        self.patient_graph[math.ceil(self.n / 2)] \
            [math.ceil(self.m / 2)] = 1

        # Initialize our counters for each patient type
        self.It = []
        self.St = []
        self.Rt = []
        self.Vt = []

    def run_multiple_sim(self, runs = 1000, prt = False, nu_0 = 0, nu_f = 1, increment = 0.1):
        '''
        Note: only tries different nu values currently

        Creates a member object named "sim_data"
            - This is a pandas dataframe that contains the results of
            the runs; see col_names below for basic structure
        '''
        
        col_names = ["nu_run", "individual_run", "nu", "I"]
        self.sim_data = pd.DataFrame(columns = col_names)
        # Initialize a pandas dataframe to store our results in

        nu_vals = np.arange(start = nu_0, stop = (nu_f + increment), step = increment)

        total_index = 0

        for i in range(0, len(nu_vals)):
            for j in range(0, runs):

                self.nu = nu_vals[i]
                self.run_simulation(prt)

                # Need a dictionary with new data for this run
                new_row = { 'nu_run': i, 
                            'individual_run': j, 
                            'nu': self.nu, 
                            'I': self.Rt[-1] }

                self.sim_data = self.sim_data.append(new_row, ignore_index = True)
                # Adding new data to the dataframe

                self.reset_internal()

                total_index += 1

    # Runs the simulation for a given number of runs
    def run_simulation(self, prt = False):
        '''
        model.run_simluation(True)

        Public method that runs the model once with the given parameters.

        prt: specifies whether to print the patient graph after each 
        iteration.
            - True: print graph
            - False: do not print graph

        '''
        self.days = self.__monte_carlo_sim(prt)
        return self.days
        
    # Plots the data
    def plot_data(self, plot_title = "Infection Model", save = False, show = True):

        '''
        model.plot_data(save, show)
        Plots the given data stored in the infection_model object.

        Usage:
        plot_title: Title to be shown on plot

        save: Option to save plot as png
        - True = save the plot
        - False = don't save the plot
        - Note: file saved with a default file name. 
            Name structure: infection_model_m<m value>_n<n value>_k<k value>.png

        show: Option to show the plot
        - True = show the plot
        - False= don't show the plot
        '''

        plt.title(plot_title)
        plt.xlabel("Days")
        plt.ylabel("Percent of Population")

        xval = range(0, self.days + 1)

        plt.plot(xval, self.It, "-", color = "red", label = "Infected")
        plt.plot(xval, self.St, "-", color = "blue", label = "Susceptible")
        plt.plot(xval, self.Rt, "-", color = "green", label = "Recovered")
        if self.vaccinate:
            plt.plot(xval, self.Vt, "-", color = "purple", label = "Vaccinated")
        plt.legend(loc = 0)

        # Code taken from GitHub - 
        # https://stackoverflow.com/questions/11373610/save-matplotlib-file-to-a-directory
        #   Saves the figure to local directory
        if (save):
            file_name = "infection_model_m{:d}_n{:d}_k{:d}.png".format(self.m, self.n, self.k)
            plt.savefig(file_name)

        if (show):
            plt.show()

    def plot_histograms_infected(self, save = False, show = True):
        '''
        Plots histograms for challenge 19-5
        '''

        # Get number of plots to make
        unique_nu = self.sim_data.nu.unique()

        fig, ind_plots = plt.subplots( math.ceil(len(unique_nu) / 2), 2 )
        
        fig.suptitle("Histograms of Infected Individuals")        

        for i in range(0, len(unique_nu) ):

            row = int(i / 2)
            col = i % 2

            rows_interest = self.sim_data.loc[ self.sim_data["nu"] == unique_nu[i] ]
            I_vals = list(rows_interest["I"])

            # Plot histograms
            n, bins, patches = ind_plots[row, col].hist(x = I_vals, bins = 100, \
                 color = "green", range = (0, 1))
                # n gets the data returned from the histogram function

            # Set x and y labels
            ind_plots[row, col].set( xlabel = "Proportion Infected", ylabel = "Frequency")


            y_spot = max(n) * 0.9 # Get location for top 
            y_decr = max(n) * 0.1 # Get decrement for text values on graph
            
            # Show nu value on plot:
            ind_plots[row, col].text(x = 0.6, y = y_spot, s = "nu = {:3f}".format(unique_nu[i]))

            # Show mean on plot:

            mean_I = np.mean(I_vals)

            ind_plots[row, col].text(x = 0.6, y = y_spot - y_decr, \
                s = "Mean = {:3f}".format(mean_I))

            ind_plots[row, col].axvline(x = mean_I, color = 'r', \
                linestyle = 'dashed')
            
            # Show variance on plot:

            var_I = np.var(I_vals)

            ind_plots[row, col].text(x = 0.6, y = y_spot - y_decr * 2, \
                s = "Var = {:3f}".format(var_I))

        if save: # Save plot if specified
            fname = "infected_histogram"
            plt.savefig(fname)

        if show: # Show plot if specified
            plt.show()

    def generate_dataframe(self):
        '''
        dataframe = model.generate_dataframe()

        Outputs a pandas dataframe containing information currenly held in object
        - Will not output dataframe if the model has not yet been ran
        '''

        len_data = [ len(self.It), len(self.St), len(self.Rt), len(self.Vt)]

        if(0 in len_data):
            print("Run model before creating dataframe")

        else:
            col_names = [ "day", "Infected", "Susecptible", "Recovered", "Vaccinated"]

            df = pd.DataFrame(columns = col_names)
            df["day"] = range(0, self.days)
            df["Infected"] = self.It
            df["Susceptible"] = self.St
            df["Recovered"] = self.Rt
            df["Vaccinated"] = self.Vt

            return df

    def current_parameters(self):
        '''
        parameters = model.current_parameters()

        Returns a dictionary containing the current parameters for the model. 
        '''
        parameters = {}

        parameters["m"] = self.m
        parameters["n"] = self.n
        parameters["k"] = self.k
        parameters["tau"] = self.tau
        
        if self.mobility:
            parameters["delta"] = self.delta
        
        if self.vaccinate:
            parameters["nu"] = self.nu

    
