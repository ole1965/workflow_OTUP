function [pulseGxyStep] = prepGzStep(kxyShift,...
    ampGradMax,rateSlewMax,tDelta,factorOs,tDuration)

if nargin < 6, tDuration = 0; end
if nargin < 5, factorOs = 1; end
tDeltaFunc = tDelta*factorOs;

if abs(kxyShift) < 1e-6, pulseGxyStep = []; return; end

if tDuration == 0
    [agxyShift,rsxyShift,nxtShift] = ...
        paraGrCpt(kxyShift,ampGradMax,rateSlewMax,tDeltaFunc); 
else
    [agxyShift,rsxyShift,nxtShift,kxyShiftRenew] = ...
        paraGrCpt(kxyShift,ampGradMax,rateSlewMax,tDeltaFunc,tDuration); 
end
pulseGxyStep = pulseGrCpt(agxyShift,rsxyShift,nxtShift,tDeltaFunc);

pulseGxyStepFunc = zeros(factorOs,length(pulseGxyStep));
for indexFactorOs = 1:factorOs
    pulseGxyStepFunc(indexFactorOs,:) = pulseGxyStep;
end
pulseGxyStep = reshape(pulseGxyStepFunc,[],1);

end