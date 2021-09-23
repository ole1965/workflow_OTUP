function [pulseRfPerCoilam,error,sta] = pulseRFcalcLsqr(mapB1,mapB0,dpDimXYZ3D,phiTarget,pulseGrTotal,tDelta,maxiter)
% Spatial Domain Method from Grissom et al.* implemented by Tingting Shao
% *= Grissom W, Cy Y, Zhang Z, Stenger VA, Fessler JA, Noll DC. Spatial domain method for the design of RF pulses in multicoil parallel excitation. Magn Reson Med. 2006;56:620-629.

    gamma = 2.6752218744e8;
    nSpatialSample = size(dpDimXYZ3D,1);
    %% RF pulse calculation
    pulseGrFunc = pulseGrTotal;
    kxyFuncSum = [cumsum(pulseGrTotal)]*gamma/2/pi*tDelta;
    % kxyFuncShift = kxyFuncSum(1) - kxyFuncSum(end);
    kxyFuncShift = kxyFuncSum(1,:) - kxyFuncSum(end,:);
    kxyFunc = kxyFuncShift + [cumsum(pulseGrFunc)]*gamma/2/pi*tDelta;
    kTrjtMinus = gamma/2/pi*tDelta/2*([pulseGrFunc]);
    kTrjt = kxyFunc - kTrjtMinus; nKxyFunc = length(kTrjt);
    % kTrjt = cumsum(pulseGrTotal)*(gamma/(2*pi)*tDelta);
    % nKxyFunc = length(kTrjt);
    %%
    clear mskTarget mzTarget  kxyFuncBevor pulseGrFunc pulseGxyFuncBevor infoGr kTrjtMinus kxyFunc kxyFuncBevor 
    if mapB0==0
        matrixA = -1i*2*pi*(real(dpDimXYZ3D(:,1)*conj(kTrjt(:,1).')) + dpDimXYZ3D(:,2)*(kTrjt(:,2).'));
    else    
    %     tIntgr = indiceKxyFunc*tDelta;
        tIntgr = ((1:nKxyFunc)-nKxyFunc)*tDelta;
        matrixA = -1i*2*pi*(real(dpDimXYZ3D(:,1)*conj(kTrjt(:,1).')) + dpDimXYZ3D(:,2)*(kTrjt(:,2).') + mapB0*tIntgr);
    end
    clear kTrjtMinus kxyFunc kxyFuncBevor  indiceKxyFunc tIntgr
    matrixA = 1i*gamma*tDelta*exp(matrixA);
    nChannelTx = size(mapB1,2); nKTrjtValid = size(kTrjt,1);
    matrixAfull = zeros(nSpatialSample,nKTrjtValid*nChannelTx);
    for indexChannelTx=1:nChannelTx
        for indexSampleTime=1:nKTrjtValid
            matrixAfull(:,indexSampleTime+nKTrjtValid*(indexChannelTx-1)) =  matrixA(:,indexSampleTime).*mapB1(:,indexChannelTx);
        end
    end
    clear matrixA kTrjt
    tol=1.e-13;
    [matrixB2,~] = lsqr(matrixAfull,phiTarget,[],maxiter);
    % tic;matrixB2 = magnitude_lsqr(matrixAfull,phiTarget,20);toc;
    sta = rad2deg(abs(matrixAfull*matrixB2));
    mxyTargetamAmp = sta - rad2deg(phiTarget);
    rmse = sqrt(sum(mxyTargetamAmp.^2)/nSpatialSample);
    error = rmse;
    %% pulse splitting
    pulseRfPerCoilam = reshape(matrixB2,size(pulseGrTotal,1),nChannelTx);
end