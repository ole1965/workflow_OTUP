function [rmse] = ZF_spins(paras,mapB1,mapB0,spatial,targetVec,maxiter)
%function to create a SPINS trajectory according the parameters from Eq. 2-4 in the manuscript.
%For the corresponding gradient waveform, the RF pulse will be calculated. 
%The output parameter 'rmse' is the rmse the RF pulse and the gradient waveform produce.
    tDelta = 1e-5;
    gamma = 2.6752218744e8;
    
    T = paras(6)*tDelta;% 0.001;
    t = 0:tDelta:T-tDelta;
    kmax = paras(1);%2;
    alpha = paras(2);%5;
    beta = paras(3);%0.5;
    u = paras(4)*2*pi;%8*pi;
    v = paras(5)*2*pi;%2*pi;
    kr = kmax./(1+exp(alpha*(t/T-beta)));
    ktheta = u*(t/T);
    kphi = v*(t/T);
    kx = kr.*cos(kphi).*sin(ktheta);
    ky = kr.*sin(kphi).*sin(ktheta);
    kz = kr.*cos(ktheta);
    k_tra = [kx(:) ky(:)  kz(:) ];

    pulseGrTotal = diff([k_tra]/(gamma/(2*pi)*tDelta));
    
    ampl = max(max(abs(pulseGrTotal)));
    sr2 = max(sqrt(sum(diff(pulseGrTotal).^2,2))/tDelta);

    if ampl > 4e-2 || sr2 > 190
            rmse = 1000;
    else
        g = [pulseGrTotal(:,1)+1i*pulseGrTotal(:,2),pulseGrTotal(:,3)];
        
        %calculate the pulse based on the current gradient waveform and put out the corresponding rmse
        [pulse,rmseSTA,~] = pulseRFcalcLsqr(mapB1,mapB0,spatial,targetVec,g,tDelta,maxiter);
        
        %Simulate with Bloch equations and calculate rmse
        mapB0 = mapB0.';
        spatial = [real(spatial(:,1)),imag(spatial(:,1)),spatial(:,2)]'; %dist = dist(:,logical(mskTarget));
        g = [real(g(:,1)),imag(g(:,1)),g(:,2)]';
        
        mag = bloch3D_mex(mapB0,spatial,pulse,mapB1,g,tDelta);
        mz = mag(3,:).';
        mzDeg = acosd(mz);
        rmse = rad2deg(targetVec) - mzDeg;
        rmse = sqrt(sum(rmse.^2)/length(targetVec));
    end
end
