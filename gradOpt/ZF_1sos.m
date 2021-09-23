function [rmse] = ZF_1sos(paras,mapB1,spatial,targetVec,maxiter,mapB0)
%function to create a single spiral according the parameters from Eq. 1 in the manuscript.
%The spiral has a specific position on the z-axis of the k-space (paras(5)).
%The start and the end of the spiral will be connected with the k-space center (according to the hardware limits). 
%For the corresponding gradient waveform, the RF pulse will be calculated. 
%The output parameter 'rmse' is the rmse the RF pulse and the gradient waveform produce.
    
    srStep = 134.35;%max. slew rate for each gradient, respectivly
    gamma = 2.6752218744e8;
    tDelta = 1e-5;
    
    %create the 2D spiral in k-space and the corresponding gradient waveform
    theta = linspace(1,0,paras(4)).';
    ga = paras(1);
    alpha = paras(2)';
    nTurns = paras(3);
    omega = 2*pi*nTurns;
    vd1 = ga*theta.^alpha.*exp(1i*omega.*theta); %Eq. 1
    v = [real(vd1),imag(vd1)];
    
   	G = diff([v]/(gamma/(2*pi)*tDelta));
    
    %check if amplitude und slew rate limits are satisfied
    ampl = max([max(G(:,1)),max(G(:,2))]);
    sr = max(sqrt(sum(diff([G]).^2,2))/tDelta);
    if ampl > 4e-2 || sr > 190
        rmse = 1000;
        return;
    end
    
    %create a path from k-space center to start of the spiral ...
    vd1 = vd1.*(cos(2*pi-angle(vd1(1)))+1i*sin(2*pi-angle(vd1(1)))); %rotate the 2D spiral so that angle(vd1(1)) = 0
    shift = vd1(1);
    vd1 = vd1 - ((vd1(1)));
    pulseGxySpiral = diff([vd1]/(gamma/(2*pi)*tDelta));
    gxyStep = pulseGxySpiral(1);
    pulseGxyStep = prepGxyStep(shift,gxyStep,0.04,srStep,tDelta,1,0);
    nGxyStep = length(pulseGxyStep);

    pulseGzStep = prepGzStep(paras(5),0.04,srStep,tDelta,1);
    nGzStep = length(pulseGzStep);
    
    if nGxyStep >= nGzStep
        Step = [pulseGxyStep,zeros(nGxyStep,1)];
        Step(end-nGzStep+1:end,2) = pulseGzStep;
        G = [Step;[pulseGxySpiral,zeros(size(pulseGxySpiral,1),1)]];
    else
        Step = [zeros(nGzStep,1),pulseGzStep];
        Step(end-nGxyStep+1:end,1) = pulseGxyStep;
        G = [Step;[pulseGxySpiral,zeros(size(pulseGxySpiral,1),1)]];
    end
    
    %... and from the end of the spiral back to center
    Gende = [[real(G(end,1)),imag(G(end,1)),G(end,2)];[0,0,0]];
    slew = sqrt(sum(diff(Gende).^2,2))/tDelta;pos = size(Gende,1)-1;
    if slew>190
        while pos>0
            sr = slew(pos);
            if sr>190
                Gende=Gende*1e3;
                difference = Gende(pos,:) - Gende(pos+1,:);
                g2 = Gende(pos,:) - difference * 1/2;% * 2/3;
                g1 = Gende(pos,:)- g2;
                Gende = [Gende(1:pos-1,:);g1;g2;Gende(pos+1:end,:)]; %insert that in the excisiting gradient shape
                Gende=Gende/1e3;
                % grad vector is longer, so we have to check from here on
                % (including the point we checked already)
                pos = pos+2;
                slew = sqrt(sum(diff(Gende).^2,2))/tDelta;
                if pos>length(slew)
                    pos = length(slew);
                end
            else
                pos = pos-1;
            end
        end
        Gende = [Gende(:,1)+1i*Gende(:,2),Gende(:,3)];
        G = [G(1:end-1,:);Gende];
    end
  
    pulseGzStep = prepGzStep(-paras(5),0.04,srStep,tDelta,1);
    nGzStep = length(pulseGzStep);
    Step = [zeros(nGzStep,1),pulseGzStep];
    G = [G;Step];
    
    %check if maximum length of 10ms is satisfied and again if slew rate limits are satisfied 
    if size(G,1)>1000 || max(sqrt(sum(diff([real(G(:,1)),imag(G(:,1)),G(:,2)]).^2,2))/tDelta) > 190
        rmse = 1000;
        return;
    end
    g = G;
    
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