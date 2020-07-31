# Infectious Disease Modeling

This repository contains code for a project that I completed for Dr. Michael Berry at UTK in COSC 377 - "Honors Introduction to Scientific Computing" during the summer of 2020. This project was based on [the writeup](https://github.com/owencqueen/infection_modeling/blob/master/Infection_sim.pdf) contained in this repo. It is important to note that this repo only contains challenges 19.1 - 19.5.

For this project, I was asked to build a simple model that would demonstrate the dynamics of an infectious disease on a network of people. This project seemed quite appropriate given the COVID-19 outbreak that was currently sweeping the US. This project explores how isolation strategies and vaccination processes can affect the extent of an outbreak. 

This project is broken down into five sub-repos, one for each section of the writeup. In each of these sub-repos, you will find a README containing descriptions of my solutions to the problem as well as my findings from the investigations. Each of the sub-repos are linked below:

  1. **[Challenge 19.1](https://github.com/owencqueen/infection_modeling/tree/master/challenge_19-1):** Model setup and examination of outbreak statistics
  2. **[Challenge 19.2](https://github.com/owencqueen/infection_modeling/tree/master/challenge_19-2):** Introduce mobility for patients
  3. **[Challenge 19.3](https://github.com/owencqueen/infection_modeling/tree/master/challenge_19-3):** Examine our model over many trials
  4. **[Challenge 19.4](https://github.com/owencqueen/infection_modeling/tree/master/challenge_19-4):** Explore how vaccination rates affect the extent of the outbreak
  5. **[Challenge 19.5](https://github.com/owencqueen/infection_modeling/tree/master/challenge_19-5):** Find optimal vaccination rates to contain the outbreak

## Running the code
Each sub-repo is setup to include only the scripts, modules required for that challenge, and the resulting plots generated from the analysis. The scripts are contained in files with the name `run_challenge_19-<challenge_number>.py` where `<challenge_number>` corresponds to the section of challenge 19 in the writeup. The module `infection_model.py` contains all of the code needed to run the scripts (other than local modules). 
