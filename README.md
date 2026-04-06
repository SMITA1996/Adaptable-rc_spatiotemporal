Adaptable Reservoir Computing Model for prediction of spatiotemporal tipping by Smita Deb and Zhengmeng Zhai.

This is the repository for our preprint titled "Anticipating tipping in spatiotemporal systems with machine learning". 
This research focuses on using an adaptable reservoir Computing (RC) model to predict spatiotemporal tipping points for various systems, including CMIP5 climate projections.

# User instructions
Compilation of the given codes requires MATLAB version R2020b and Fortran.

# Codes and results
adaptable_rc_spatiotemporal.m can be run to obtain predictions of adaptable RC on NMF-reduced data.

The folder data generation contains two subfolders that contain code for generating spatial patterns for the continuous model and the discrete cellular automata model.

numerical_data.f90 (contained in the subfolder VT) can simulate the sequence of spatial patterns along the bifurcation gradient for all other parameters kept fixed for the vegetation turbidity model.
The  code can be modified to obtain data for other continuous models by changing the deterministic equation and corresponding parameters.

The code ca.m in the subfolder ca can simulate the sequence of spatial patterns along the bifurcation gradient for the discrete (cellular automata) model.





