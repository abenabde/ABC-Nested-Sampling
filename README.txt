%-------------------------------------------------------------------------------------------------------------------------------------------------------
% Example 1 : SDOF system
% ABC-NS for  parameter estimation
%-------------------------------------------------------------------------------------------------------------------------------------------------------
The single-degree-of-freedom (SDOF) system governed by the following ordinary differential equation is considered : 

                                                              mx'' + cx' + kx =f(t) 

where m, c and k are the mass, the damping and the stiffness parameters 
The mass is supposed to be known 
f(t) : is the excitation force
The Normalised Mean Square Error (NMSE) is used to measure the discrepancy between the training data  and the 
 simulated data. 

For more details about the implementation and the choice of the tuning parameters, please refer to the references [1,2] and 
 the slides:
https://figshare.shef.ac.uk/articles/Recent_Advances_in_Approximate_Bayesian_Computation_Methodology_Application_in_structural_dynamics_/4578901

To get an idea about the principle of the Nested Sampling algorithm, please refer to the references [4,5].
%-------------------------------------------------------------------------------------------------------------------------------------------------------
% The folder 'DEMO_PE' contains the following scripts: 
%-------------------------------------------------------------------------------------------------------------------------------------------------------
- Run.m: run the ABC-NS algorithm for parameter estimation
- model.m: the script simulate a data from the model for a given set of parameters
- force.m: the external force 
- fun.m: solve the ordinary differential equation
- draw_from_ellipsoid.m: draws points uniformly from an ndims-dimensional ellipsoid
- calc_ellipsoid.m : calculate properties of ellipsoid given a set of particles

Notice : 
The scripts  'calc_ellipsoid.m' and 'draw_from_ellipsoid.m'  are taken from the MultiNest package [3]. 
%-------------------------------------------------------------------------------------------------------------------------------------------------------
% References: 
%-------------------------------------------------------------------------------------------------------------------------------------------------------
[1] A. Ben Abdessalem, N. Dervilis, D. Wagg, K. Worden, Model selection
and parameter estimation of dynamical systems using a novel variant of approximate
Bayesian computation, Mechanical Systems and Signal Processing. 122 (2019) 364-386.

[2] A. Ben Abdessalem, N. Dervilis, D. Wagg, K. Worden, An Efficient Likelihood-Free Bayesian Computation
 for Model Selection and Parameter Estimation Applied to Structural Dynamics. 
Structural Health Monitoring, Photogrammetry & DIC, Volume 6 pp 141-151.

[3] F. Feroz, M.P. Hobson and M. Bridges. MULTINEST: an efficient and robust Bayesian inference tool for cosmology and particle
physics. Mon. Not. R. Astron. Soc., 398:1601-1614, 2009.

[4] J. Skilling (2006) Nested sampling for general Bayesian computation.
 Bayesian Analysis, 1(4):833-860 [http://dx.doi.org/10.1214/06-BA127].

[5] P. Mukherjee, D. Parkinson and A.R. Liddle (2006) A Nested Sampling Algorithm for 
Cosmological Model Selection. Astrophysical Journal Letters, 638 (2). L51-L54 [http://dx.doi.org/10.1086/501068]
%-------------------------------------------------------------------------------------------------------------------------------------------------------
% FeedBack 
%------------------------------------------------------------------------------------------------------------------------------------------------------- 
If you find any bug or have any suggestion, please send an email to :

                                      anis.ben-abdessalem@univ-angers.fr
%-------------------------------------------------------------------------------------------------------------------------------------------------------