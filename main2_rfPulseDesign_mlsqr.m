%Script to perform the UP design with the magnitude least-squares optimization as described in the manuscript and step 2 in the 'readMe.md'.
setUp;
%create design database
mapB1 = []; mapB0 = []; mask = []; spatial = [];
disp('Prepare design database...')
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

%Select underlying transmit gradient waveform
disp('Choose underlying gradient waveform for RF pulse design...');
[gradFile,path] = uigetfile('gradientWaveforms\*.mat','Choose underlying gradient waveform for RF pulse design.');
load([path,gradFile]);

disp('Design RF Pulse...');
maxiter = 20;
%Least-squares optimization (spatial domain method)
[pulseRf,err,sta] = pulseRFcalcLsqr(mapB1,mapB0,spatial,targetVec,gr,tDelta,maxiter);
% %Magnitude-Least-squares optimization
% [pulseRf,err,sta] = pulseRFcalcMLsqr(mapB1,mapB0,spatial,targetVec,gr,tDelta,maxiter);

disp(['Amplitude of resulting pulse: ',num2str(max(abs(pulseRf(:)))),' V'])

lVec = sizeMatrix(1)*sizeMatrix(2)*sizeMatrix(3);

gr = [real(gr(:,1)),imag(gr(:,1)),gr(:,2)]';
disp('Simulate Pulse on test database...');
for j = 1:length(headsTest)
    disp(['Simulate on head ',num2str(headsTest(j))])
    switch target
        case 1 %wB
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_wB(['B0B1datasets\recon_subj',num2str(headsTest(j)),'.mat']);
        case 2 %LEx
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj',num2str(headsTest(j)),'.mat']);
            %postprocess the design mask according to the nuclei application in the paper
            masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
            masksingle(~mask2) = false;
        case 3 %LEx
            [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj',num2str(headsTest(j)),'.mat']);   
            masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
    end 
    masksingle = logical(masksingle(:));
    %Prepare the Data for the Mex Bloch Simulation
    mapB0single = mapB0single(logical(masksingle))';
    spatial = [real(spatialsingle(logical(masksingle),1)),imag(spatialsingle(logical(masksingle),1)),spatialsingle(logical(masksingle),2)]';
    
    mapB1single = mapB1single(logical(masksingle),:);

    mag = bloch3D_mex(mapB0single,spatial,pulseRf,mapB1single,gr,tDelta);    
    mz = mag(3,:).';
    mzDeg = acosd(mz);
    
    targetVecSingleMasked = targetVecSingle; targetVecSingleMasked(~masksingle) = [];
    rmse = rad2deg(targetVecSingleMasked) - mzDeg;
    rmse = sqrt(sum(rmse.^2)/sum(masksingle));
    nrmse = rmse/fa;
    
    mzplot = zeros(lVec,1);
    mzplot(logical(masksingle)) = mzDeg;
    
    slTra = 36;slSag = 31;slCor = 31;
    M = (reshape(mzplot,sizeMatrix));
    M1 = M(:,:,slTra);M2 = rot90(squeeze(M(:,slSag,:)));M3 = rot90(squeeze(M(slCor,:,:)));
    MM = [rot90(M1);M2;M3];
    figure;b = imagesc(MM);set(b,'AlphaData',MM~=0);set(gca,'xtick',[]);set(gca,'ytick',[]);hold on;
    axis equal tight;caxis([0 fa]);
    title(['Head ',num2str(headsTest(j)),', RMSE: ',num2str(rmse),' , NRMSE: ',num2str(nrmse)]);
end

%save the resulting pulse in folder rfPulses
savePulseRf;

