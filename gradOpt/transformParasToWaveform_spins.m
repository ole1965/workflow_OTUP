function gr = transformParasToWaveform_spins(paras)
    tDelta = 1e-5;
    T = paras(6)*tDelta;
    t = 0:tDelta:T-tDelta;
    kmax = paras(1);
    alpha = paras(2);
    beta = paras(3);
    u = paras(4)*2*pi;
    v = paras(5)*2*pi;
    kr = kmax./(1+exp(alpha*(t/T-beta)));
    ktheta = u*(t/T);
    kphi = v*(t/T);
    kx = kr.*cos(kphi).*sin(ktheta);
    ky = kr.*sin(kphi).*sin(ktheta);
    kz = kr.*cos(ktheta);
    k_tra = [kx(:) ky(:)  kz(:) ];
    
    gamma = 2.6752218744e8;
    pulseGrTotal = diff([k_tra]/(gamma/(2*pi)*tDelta));

    gr = [pulseGrTotal(:,1)+1i*pulseGrTotal(:,2),pulseGrTotal(:,3)];
end