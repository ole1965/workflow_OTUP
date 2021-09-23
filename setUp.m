clear;
addpath('general');
addpath('blochMex');
addpath('gradientWaveforms');

%% Define the target
%Choose the target excitation pattern and the desired flip angle. The prepared patterns are the whole brain excitaion pattern
%and the nuclei local excitation pattern from the manuscript.
fa = 7;
target = 2; %1 for whole brain (wB) Ex,
            %2 for LEx of targetNuclei (according to the nuclei exc. patern from the paper),
            %3 for LEx of targetM (according to the M exc. patern from the paper), 
%% Define the slices of interest
%Set the transversal slices which shall be part of the design mask (only important for the local excitation application in the manuscript)
%The transversal slices from 'startPos' to 'endPos' will be included in the design mask, the other slices will not be included.
%targetWB: startPos = 29;endPos = 44;
%targetNuclei: startPos = 29;endPos = 36;
%targetM: startPos = 36;endPos = 36;
startPos = 29;endPos = 36;    
nrSl = endPos-startPos+1; %number of included slices

%% Define database and test heads
%Choose the heads to put inside the pulse design database and the heads on which the pulse should be tested.
headsDB = [1:3];
headsTest = [1];

%% Prepare the target matrix and vector
[targetMat,mask2] = prepare_target(target,startPos,endPos,nrSl);
targetMat(targetMat>0) = deg2rad(fa);
targetVec = targetMat(:);



%% Set up miscellaneous variables
sizeMatrix = [64,64,64];
lVec = sizeMatrix(1)*sizeMatrix(2)*sizeMatrix(3);
name = ['DB',num2str(length(headsDB)),'_FA',num2str(fa),'_'];
tDelta = 1e-5;
gamma = 2.6752218744e8;

%% settings for databaseSize Investigations
nrDatasets = 3;