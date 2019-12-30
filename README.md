# MCSS2020_HyperbolicISS
Code for the examples presented in a paper submitted at MCSS
The repository contains two folders, one for each of the two examples presented in the paper. 
Folder Example2 contains two subfolders, i.e., Example2.1 and Example2.2, each containing the
two scenarios considered in Example 2. 

For each of the examples, the code needs to be used as follows: 

1- Run the file Design.m, this defines the plant data and designs the controller gain via the results presented in the paper.  
2- Go into the subfolder Simulaiton and run the file Run.m. This simulates the closed-loop system and generates the components of plant trajectoty x(t,z) in the form of a matrix compatible with "surf" and "mesh" instructions.     
