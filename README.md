# Closed-form-Power-Flow

This file contain the code related to the paper: 
Parikshit Pareek and Hung Nguyen, "A Framework for Analytical Power Flow Solution using Gaussian Process Learning" IEEE Transactions on Sustainable Energy, Sept. 2021. 
Link: https://ieeexplore.ieee.org/document/9552521
Preprint:  https://www.researchgate.net/publication/354937182_A_Framework_for_Analytical_Power_Flow_Solution_using_Gaussian_Process_Learning

Cite As: 
```
@article{pareek2021framework,
  title={A Framework for Analytical Power Flow Solution using Gaussian Process Learning},
  author={Pareek, Parikshit and Nguyen, Hung D},
  journal={IEEE Transactions on Sustainable Energy},
  year={2021},
  publisher={IEEE}
}
```

In perticular, the code can be used to obtain the closed-form power flow expression and performance analysis results as given in table 1 and section V-A of the manuscript. 

## Details of Files: 
`CFPF_multikernel.m` : Main file to run the code for obtaining the CFPF Approximation and Performance for different kernels.
`input_dataset_Load.m` : Creating load data set for training and testing 
`MCS_output.m`    : Monte-Carlo Simulation to obtain testing data points
`rand_sample_x.m` : Generating random samples
`runpf_complete.m` : MATPOWER codes combined together to avid dependencies 
`Sampling_Jaco.m`  : Cover to 'runpf' for obtaining power flow datasets

Dependencies: 
GAUSSIAN PROCESS REGRESSION AND CLASSIFICATION Toolbox version 4.2 for GNU Octave 3.2.x and Matlab 7.x and higher.
Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-06-11.
Link: http://www.gaussianprocess.org/gpml/code/matlab/doc/

Newer GPML toolbox might require matching hyperparmeter initialization.

Variables and Code Outputs: 

1. Loutkr : Cell array containing structre of results for each kernel

N_train: Number of training samples
N_bus : Number of PQ buses
N_test: Number of testing samples 
#Hyper_parameters : Number of hyper-parameters
D   : Number of random injections

     alphaV: [N_train × N_bus double] : Value of alpha for Voltage Magnitude from equation (4)
        muV: [N_test × N_bus double]  : Mean prediction of Voltage Magnitude from equation (5)
        s2V: [N_test × N_bus double]  : Predicitive variance of Voltage Magnitude from equation (16)
      sf_lV: [#Hyper_parameters × N_bus double] : Hyper-parameters for Voltage Maginitide learning
        ytV: [N_train × N_bus double]: Training samples of Voltage Maginitide

    alphaTh: [N_train × N_bus double] : Value of alpha for Voltage Angle from equation (4)
       muTh: [N_test × N_bus double]   : Mean prediction of Voltage Angle from equation (5)
       s2Th: [N_test × N_bus double]   :  Predicitive variance of Voltage Angle from equation (16)
     sf_lTh: [#Hyper_parameters × N_bus double] : Hyper-parameters for Voltage Angle learning
       ytTh: [N_train × N_bus double] :  Training samples of Voltage Angle

2. maeV (maeTh) : Mean Absolute Error Voltage Magnitude (Angle)

3. Mcskr : Cell array for MCS Results


          ` V: [N_test × N_bus double double]
      erV_par: [N_test × N_bus double double]
       erV_L1: 0.037771 % L_1 Norm Error in |V|
       erV_L2: 0.024608 % L_2 Norm Error in |V|
     erV_Linf: 0.062861 % L_inf Norm Error in |V|`

         Thac: [N_test × N_bus double double]
     erTh_par: [N_test × N_bus double double]
      erTh_L1: 2.6544 % L_1 Norm Error in Voltage Angle 
      erTh_L2: 3.5956 % L_2 Norm Error in Voltage Angle 
    erTh_Linf: 7.7006 % L_inf Norm Error in Voltage Angle 



4. D_learningkr: Learning data cell array, each array has a structre with data
                    xx: [N_train × D double] Training dataset or Design Matrix
                    xs: [N_test × D double]  Testing dataset 
Other variable names are self explaintory
