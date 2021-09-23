function [pulseG] = pulseGrCpt(ampGrad,rateSlew,nSampleTime,tDelta)
% ampGrad=4;rateSlew=15;tDelta;nSampleTime=32*paraTbwScale+1;
% load aa.mat;
% gamma = 2*pi*4.257;     % krad/G

if nSampleTime*ampGrad == 0, pulseG = []; return; end

if imag(ampGrad)==0 || real(ampGrad)==0
    pulseG = zeros(nSampleTime,1);          % nSample to be 2*N+1
    nSampleT1 = round(ampGrad/rateSlew/tDelta); 
    nSampleT2 = round(nSampleTime-1-nSampleT1*2);
    pulseG(1:nSampleT1) = [0:nSampleT1-1]*rateSlew*tDelta;
    pulseG(nSampleT1+1:nSampleT1+nSampleT2+1) = ampGrad;
    pulseG(nSampleT1+nSampleT2+2:nSampleTime) = ...
        [nSampleT1-1:-1:0]*rateSlew*tDelta;
%     kTrjt = gamma/2/pi*cumsum(pulseG)*tDelta;
%     figure;plot((0:length(pulseG)-1)*tDelta,pulseG); figure;plot(kTrjt);
else
    ampGradx = real(ampGrad); rateSlewx = real(rateSlew);
    pulseGx = zeros(nSampleTime,1);          
    nSampleT1 = round(ampGradx/rateSlewx/tDelta); 
    nSampleT2 = round(nSampleTime-1-nSampleT1*2);          
    pulseGx(1:nSampleT1) = [0:nSampleT1-1]*rateSlewx*tDelta;
    pulseGx(nSampleT1+1:nSampleT1+nSampleT2+1) = ampGradx;
    pulseGx(nSampleT1+nSampleT2+2:nSampleTime) = ...
        [nSampleT1-1:-1:0]*rateSlewx*tDelta;
    ampGrady = imag(ampGrad); rateSlewy = imag(rateSlew);
    pulseGy = zeros(nSampleTime,1);          
    nSampleT1 = round(ampGrady/rateSlewy/tDelta); 
    nSampleT2 = round(nSampleTime-1-nSampleT1*2);
    pulseGy(1:nSampleT1) = [0:nSampleT1-1]*rateSlewy*tDelta;
    pulseGy(nSampleT1+1:nSampleT1+nSampleT2+1) = ampGrady;
    pulseGy(nSampleT1+nSampleT2+2:nSampleTime) = ...
        [nSampleT1-1:-1:0]*rateSlewy*tDelta;
    pulseG = pulseGx + 1i*pulseGy;   
    
%     figure;subplot(2,1,1);plot((0:length(pulseG)-1)*tDelta,pulseGx);
%     subplot(2,1,2);plot((0:length(pulseG)-1)*tDelta,pulseGy); 
%     kTrjt = gamma/2/pi*cumsum(pulseG)*tDelta;
%     figure;plot(real(kTrjt),imag(kTrjt));
end

% disp(kTrjt(end));

end