# workflow_gradOptimization_UPdesign

Titel of Paper 3

As described in the manuscript (title), in this scripts the following calculation steps are implemented: transmit gradient waveform optimization and
UP design using magnitude-least-square optimization or/and UP design fmincon and bloch equations. Each calculation step can be run independently from each other and
is implemented in an own script.
At the begin of each script the 'setUp.m' script is executed. Herein the main settings are defined. These main settings are (among others):
- desired target excitation pattern
- desired target flip angle
- desired transversal slices of interest to include in the design mask
- set which heads to include in the design database and on which heads the pulse should be tested (important for UP calculation).
The three example target patterns from the manuscript (targetWB, targetNuclei and targetM) are provided.
The B1, B0 and further datasets from three example heads can be downloaded from: https://www.dropbox.com/sh/6jlfdk5cz2ej7fo/AADAXuxOT5jA3L7Z32M4lcqXa?dl=0 
New datasets and target excitation pattern can be easily added.

Detailed description of the four calculation steps:

1. 	Optimization of the transmit gradient waveform
    four different basis trajectories were chosen:
	A single variable density (vd) spiral-in trajectory (1SOS), a two stack of vd spiral-in trajectory (2SOS),
	a three stack of vd spiral-in trajectory (3SOS) and a spiral nonselective trajectory (SPINS).
	For each of the four basis trajectories there is a extra optimization script: 'main1_gradOpt_1sos.m', 'main1_gradOpt_2sos.m', 'main1_gradOpt_3sos.m' and 'main1_gradOpt_spins.m'.
	In each script the corresponding basis trajectory is optimized using the particle swarm optimization.
    These calculations can be done on a pc or cluster (recommanded).
	The corresponding submit scripts for a slurm cluster can be found in the folder 'slurmSubmitBatchScripts'.
	The optimized parameter sets for each basis trajectorie are saved at the end of each script (folder: 'gradientWaveforms\optResults').
	The user can check the resulting rmse values of each parameter set in the 'rmse_val' variable. 
	The gradient waveform that provides the minimum rmse is saved in the folder 'gradientWaveforms'.
	These scripts can be used to optimize gradients not only for UP applications but also for TP design. 
    The optimized gradient waveforms and trajectories are stored in 'gradientWaveforms\gradientWaveformsPaper'.
    
	
2.	UP calculation using the magnitude-least-square optimization.
	The script 'main2_rfPulse_Design_mlsqr.m' executes the UP rf pulse design based on a selected transmit gradient waveform (preferably an optimized one) and a design database.
	The design algorithm is the magnitude-least-square optimization. However, the least-square (spatial domain method) optimization is also provided.
	The resulting rf pulse will be tested on the heads in the test database (as set in the 'setUp.m' script) using bloch equations.
	The resulting pulse will be saved as the pulse itself and in a form to use it as initial guess in step 4 (in folder: 'rfPulses').
	In general, for small tip angles these method produces 'close to optimum' pulses. 
	
3. 	UP calculation using fmincon and the bloch equations.
	The script 'main_rfPulseDesign_fminconBloch_pc.m' executes the UP rf pulse design based on a selected transmit gradient waveform (preferably an optimized one),
    a selected start pulse and a design database. The pulse will be designed using the 'fmincon.m' function exploiting Blochequations implemented in a Mex file.
    An start pulse of initial guess is needed to start this optimization (preferably the results from 2.). 
	An output function saves the pulse and further information after each iteration. The data is stored in 'rfPulses\currentPulses_pc'.
    There is also script provided to run on a cluster ('main_rfPulseDesign_fminconBloch_cluster.m').



Miscellaneous:
	Analysis of database size influence for UP pulse performance (in the supporting information of the paper).
	The script 'main_databaseSizeInvestigation.m' creates 500 random orders for the provided B0 and B1 datasets. For each order and databasesize a rf pulse will be 
	designed using the magnitude-least-square optimization. The resulting pulses will be tested (using bloch eqautions) on a non-database head (the last one in the order).
	The resulting rmse will be stored in 'databaseSizeInvestigation\results'.
	It is recommanded to run this script in [parallel on a cluster (if possible).
	The 'evaluation.m' script to evaluate the results and create the diagramms from frigure 4 in the manuscript can be found in 'databaseSizeInvestigation'. 
