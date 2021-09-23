function [pulseGxyStep] = prepGxyStep(kxyShift,gxyShift,...
    ampGradMax,rateSlewMax,tDelta,factorOs,flagShiftInv)
% kxyShift,gxyShift,ampGradMax,rateSlewMax,tDelta

if nargin < 7, flagShiftInv = 0; end
if nargin < 6, factorOs = 1; end
tDeltaFunc = tDelta*factorOs;

if abs(kxyShift) < 1e-6, pulseGxyStep = []; return; end

gamma = 2.6752218744e8;
% 
% ampGradMax = 4;         % max gradient amplitude in G/cm
% rateSlewMax = 15;       % max slew rate in G/cm/ms, here 15 = SR150
% tDelta = 0.01;         % largest step for time sampling (2us recommended) 
% 
% kxyShift = 0.4274; nxyTurns = 5; rateSlewMaxSpiral = rateSlewMax;
% kxySpiral = trjtSpiralCpt(kxyShift,nxyTurns,ampGradMax,rateSlewMaxSpiral,tDelta); 
% kxySpiral = kxySpiral-kxySpiral(end);grSpiral = [0; ktog(kxySpiral,tDelta)];
% kxyShift = kxySpiral(1); gxyShift = grSpiral(2);

kxShift = real(kxyShift);   % let kyShift=0
[agxShift,rsxShift,nxShift]=paraGrCpt(kxShift,ampGradMax,rateSlewMax,tDeltaFunc); 
pulseGxStep = pulseGrCpt(agxShift,rsxShift,nxShift,tDeltaFunc);

gyShift = imag(gxyShift); % let gxShift=0
absGyShift = abs(gyShift); flagGyShift = gyShift/absGyShift;

nSampleTimeShiftAfter = ceil(absGyShift/rateSlewMax/tDeltaFunc)*2-1;
rateSlewShiftSr = absGyShift/nSampleTimeShiftAfter/tDeltaFunc;
pulseGyShiftSrAfter = [0:nSampleTimeShiftAfter]*rateSlewShiftSr*tDeltaFunc;
kyShiftSr = gamma/2/pi*sum(pulseGyShiftSrAfter)*tDeltaFunc;
[ampGrad,rateSlew,nSampleTimeShift] = paraGrCpt(kyShiftSr,absGyShift,rateSlewMax,tDeltaFunc);
pulseGyShiftSrAfore = pulseGrCpt(ampGrad,rateSlew,nSampleTimeShift,tDeltaFunc);
pulseGyShiftSr = flagGyShift*[-pulseGyShiftSrAfore(1:end-1); pulseGyShiftSrAfter(1:end).'];
% changed 2013.11.08
% pulseGyShiftSr = flagGyShift*[-pulseGyShiftSrAfore(1:end-1); pulseGyShiftSrAfter(2:end).'];
% changed 2013.11.08
% kyShiftSr = gamma/2/pi*cumsum(pulseGyShiftSr)*tDelta;
% figure;plot(kyShiftSr);

nyShift = length(pulseGyShiftSr);
nxyShift = max(nxShift,nyShift);
pulseGxyStep = zeros(nxyShift,1);
nxShiftInv = flagShiftInv*(nxyShift-nxShift);
nyShiftInv = flagShiftInv*(nxyShift-nyShift);
pulseGxyStep(1+nxShiftInv:nxShift+nxShiftInv) = ...
    pulseGxyStep(1:nxShift) + pulseGxStep(:);
pulseGxyStep(1+nyShiftInv:nyShift+nyShiftInv) = ...
    pulseGxyStep(1:nyShift) + 1i*pulseGyShiftSr(:);

pulseGxyStepFunc = zeros(factorOs,length(pulseGxyStep));
for indexFactorOs = 1:factorOs
    pulseGxyStepFunc(indexFactorOs,:) = pulseGxyStep;
end
pulseGxyStep = reshape(pulseGxyStepFunc,[],1);

% kxyTrjt = gamma/2/pi*cumsum(pulseGxyStep)*tDelta;

% figure;plot(real(kxyTrjt),imag(kxyTrjt));
% figure;subplot(1,2,1); plot(real(pulseGxyStep));
% subplot(1,2,2); plot(imag(pulseGxyStep));
if flagShiftInv == 1, pulseGxyStep = flipud(pulseGxyStep); end

end