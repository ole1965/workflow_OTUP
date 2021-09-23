function gr = transformParasToWaveform_2sos(paras)
    srStep = 134.35;%max. slew rate for each gradient, respectivly
    
    %spiral 1
    gamma = 2.6752218744e8;
    theta = linspace(1,0,paras(4)).';
    ga = paras(1);
    alpha = paras(2)';
    nTurns = paras(3);
    omega = 2*pi*nTurns;
    vd1 = ga*theta.^alpha.*exp(1i*omega.*theta);%Eq. 1
    
    %spiral 2
    gamma = 2.6752218744e8;
    theta = linspace(1,0,paras(8)).';
    ga = paras(5);
    alpha = paras(6)';
    nTurns = paras(7);
    omega = 2*pi*nTurns;
    vd2 = ga*theta.^alpha.*exp(1i*omega.*theta);%Eq. 1
    
    %create a path from k-space center to start of the spiral 1 ...
    vd1 = vd1.*(cos(2*pi-angle(vd1(1)))+1i*sin(2*pi-angle(vd1(1)))); %rotate the 2D spiral so that angle(vd1(1)) = 0
    shift = vd1(1);
    vd1 = vd1 - ((vd1(1)));
    pulseGxySpiral = diff([vd1]/(gamma/(2*pi)*tDelta));
    gxyStep = pulseGxySpiral(1);
    pulseGxyStep = prepGxyStep(shift,gxyStep,0.04,srStep,tDelta,1,0);
    nGxyStep = length(pulseGxyStep);
    
    pulseGzStep = prepGzStep(paras(9),0.04,srStep,tDelta,1);
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
    
    
    %... and from the end of the spiral 1 to start of spiral 2 ... 
    vd2 = vd2.*(cos(2*pi-angle(vd2(1)))+1i*sin(2*pi-angle(vd2(1)))); %rotate the 2D spiral so that angle(vd1(1)) = 0
    shift = vd2(1);
    vd2 = vd2 - ((vd2(1)));
    pulseGxySpiral = diff([vd2]/(gamma/(2*pi)*tDelta));
    gxyStep = pulseGxySpiral(1);
    pulseGxyStep = prepGxyStep(shift,gxyStep,0.04,srStep,tDelta,1,0);
    nGxyStep = length(pulseGxyStep);
    
    pulseGzStep = prepGzStep(paras(10)-paras(9),0.04,srStep,tDelta,1);
    nGzStep = length(pulseGzStep);
    
    if nGxyStep >= nGzStep
        Step = [pulseGxyStep,zeros(nGxyStep,1)];
        Step(end-nGzStep+1:end,2) = pulseGzStep;
        Gneu = [Step;[pulseGxySpiral,zeros(size(pulseGxySpiral,1),1)]];
    else
        Step = [zeros(nGzStep,1),pulseGzStep];
        Step(end-nGxyStep+1:end,1) = pulseGxyStep;
        Gneu = [Step;[pulseGxySpiral,zeros(size(pulseGxySpiral,1),1)]];
    end
    G = [G;Gneu];
    
    %... and from the end of the spiral 2 back to the center. 
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
    
    pulseGzStep = prepGzStep(-paras(10),0.04,srStep,tDelta,1);
    nGzStep = length(pulseGzStep);
    Step = [zeros(nGzStep,1),pulseGzStep];
    
    gr = [G;Step];
end