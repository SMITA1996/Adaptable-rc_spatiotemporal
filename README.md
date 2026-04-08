Adaptable Reservoir Computing Model for prediction of spatiotemporal tipping by Smita Deb and Zhengmeng Zhai.

This is the repository for our preprint titled "Anticipating tipping in spatiotemporal systems with machine learning". 
This research focuses on using an adaptable reservoir Computing (RC) model to predict spatiotemporal tipping points for various systems, including CMIP5 climate projections.

# User instructions
Compilation of the given codes requires MATLAB version R2020b and Fortran.

# Codes and results
The script adaptable_rc_spatiotemporal.m is executed to obtain predictions using adaptable reservoir computing on NMF-reduced data.


The data_generation folder contains two subfolders with code for generating spatial patterns for:
1. continuous models, and
2. Discrete cellular automata models.

Within the VT subfolder, the file numerical_data.f90 simulates sequences of spatial patterns along a bifurcation gradient while keeping other parameters fixed for the vegetation–turbidity model.
This code can be adapted to other continuous models by modifying the deterministic equations and the associated parameters.

This code is implemented in Fortran to achieve faster performance for multiple realizations. It can also be equivalently written in MATLAB or Python.

In the ca subfolder, the script ca.m generates sequences of spatial patterns along the bifurcation gradient for the discrete (cellular automata) model.










