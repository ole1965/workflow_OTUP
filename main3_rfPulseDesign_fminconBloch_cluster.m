%Script to perform the UP design with the fmincon and the bloch equations as described in the manuscript and step 3 in the 'readMe.md'.
setUp;
addpath('rfPulseDesign_fminconBloch\')

%Define underlying transmit gradient waveform
gradFile = 'gradientWaveforms\gradientWaveformsPaper\gradient_2sos_targetNuclei.mat';

%Define the start pulse for fmincon (the selected pulse need to have the same length as the selected gradient waveform)
pulseFile = 'rfPulses\pulseRf_2021-09-15T184416754DB1_FA7_2sos_targetNuclei.mat'; %that pulse may also be a resulting pulse of a earlier fmincon optimization (stored in rfPulses\currentPulseCluster)

%create design database
mapB1 = []; mapB0 = []; mask = []; spatial = [];
for j = 1:length(headsDB)
    switch target
        case 1 %targetWB
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_wB(['B0B1datasets\recon_subj',num2str(headsDB(j)),'.mat']);
            if startPos > 1 && endPos < 64
                masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
            end
        case 2 %targetNuclei
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj',num2str(headsDB(j)),'.mat']);
            %postprocess the design mask according to the nuclei application in the paper
            if startPos > 1 && endPos < 64
                masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
            end
            masksingle(~mask2) = false;
        case 3 %targetM
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj',num2str(headsDB(j)),'.mat']); 
            if startPos > 1 && endPos < 64
                masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
            end
    end
    masksingle = logical(masksingle(:));
    mapB1 = [mapB1;mapB1single];clear mapB1single;
    mapB0 = [mapB0;mapB0single];clear mapB0single;
    mask = [mask;masksingle];clear masksingle;
    spatial = [spatial;spatialsingle];clear spatialsingle
end

targetVecSingle = targetVec;
targetVec = repmat(targetVec,length(headsDB),1);

mapB1 = mapB1(logical(mask),:);
mapB0 = mapB0(logical(mask),:);
spatial = spatial(logical(mask),:);
targetVec = targetVec(logical(mask),:);
targetVecDeg = rad2deg(targetVec);

%Prepare the Data for the Mex Bloch Simulation
mapB0 = mapB0.';
spatial = [real(spatial(:,1)),imag(spatial(:,1)),spatial(:,2)]';

%Select underlying transmit gradient waveform
disp('Select underlying gradient waveform for RF pulse design...');
load([gradFile]);
gr = [real(gr(:,1)),imag(gr(:,1)),gr(:,2)]';

%Prepare the start pulse for fmincon (the selected pulse need to have the same length as the selected gradient waveform)
load([pulseFile]);
startPulse = [abs(pulseRf(:));angle(pulseRf(:))]; %Conversion into amplitude and phase (easier to constrain the pulse amplitude)

%Check if gradient waveform and rf pulse have the same length
if size(gr,2) ~= size(pulseRf,1)
    error('selected gradient waveform and start pulse have unequal durations.');
end

% Optimization
algoropt = 'active-set';
algor = 'as';

if exist('startPulse','var')
    lb = [zeros(size(startPulse(1:end/2)));-Inf*ones(size(startPulse(end/2+1:end)))];
    ub = [130*ones(size(startPulse(1:end/2)));Inf*ones(size(startPulse(end/2+1:end)))];

%     parpool([13,96]);

    %call an output function in order to save the current pulse and further information after each iteration
    options = optimset('Display','iter','Algorithm',algoropt,'UseParallel',true,'OutputFcn',str2func(['outfun_cluster']));
    disp('Starting fmincon optimization...');
    [optpuls,optval,exitflag,output] = fmincon(@(pulseRf) ZF_bloch_Mex(pulseRf,mapB1,mapB0,spatial,targetVecDeg,gr,tDelta),startPulse,[],[],[],[],lb,ub,[],options);
end