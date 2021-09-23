function [ampGrad,rateSlew,nSampleTime,kb] = paraGrCpt(ka,ampGradMax,rateSlewMax,tDelta,tDuration)
% ka=8+9*i;ampGradMax=4;rateSlewMax=15;tDelta=0.0238;

% add tDuration, kb on 20160217, nur for the imag(ka)*real(ka)==0 case

gamma = 2.6752218744e8;

if nargin < 5, tDuration = 0; end
if nargout > 3, kb = ka; end
if ka==0,ampGrad=0;rateSlew=0;nSampleTime=0;return; end

if tDuration <= 0
    
if imag(ka)*real(ka)==0
    flagKa = ka/abs(ka); ka = abs(ka);
    if ampGradMax^2/rateSlewMax >= ka/(gamma/2/pi)
        t1 = sqrt(ka/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1/tDelta); t1 = nSampleT1*tDelta;
        ampGrad = ka/(gamma/2/pi)/t1; rateSlew = ampGrad/t1; 
        nSampleTime = nSampleT1*2 + 1;
    else
        t1 = ampGradMax/rateSlewMax; t2 = ka/(gamma/2/pi)/ampGradMax-t1;
        nSampleT1 = ceil(t1/tDelta); t1 = nSampleT1*tDelta;
        nSampleT2 = ceil(t2/tDelta); t2 = nSampleT2*tDelta;
        ampGrad = ka/(gamma/2/pi*(t1+t2)); rateSlew = ampGrad/t1; 
        nSampleTime = nSampleT1*2 + nSampleT2 + 1;
    end
    ampGrad = ampGrad*flagKa; rateSlew = rateSlew*flagKa;
else    
    kx = real(ka); ky = imag(ka); caseID = 0;       
    flagKx = kx/abs(kx); kx = abs(kx); flagKy = 1i*ky/abs(ky); ky = abs(ky);
    if ampGradMax^2/rateSlewMax >= kx/(gamma/2/pi)
        t1x = sqrt(kx/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1x/tDelta); nSampleTimeX = nSampleT1*2 + 1;
        caseID = caseID + 0; 
    else
        t1x = ampGradMax/rateSlewMax; t2 = kx/(gamma/2/pi)/ampGradMax-t1x;
        nSampleT1 = ceil(t1x/tDelta); nSampleT2 = ceil(t2/tDelta); 
        nSampleTimeX = nSampleT1*2 + nSampleT2 + 1; caseID = caseID + 1; 
    end 
    if ampGradMax^2/rateSlewMax >= ky/(gamma/2/pi)
        t1y = sqrt(ky/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1y/tDelta); nSampleTimeY = nSampleT1*2 + 1;
        caseID = caseID + 2; 
    else
        t1y = ampGradMax/rateSlewMax; t2 = ky/(gamma/2/pi)/ampGradMax-t1y;
        nSampleT1 = ceil(t1y/tDelta); nSampleT2 = ceil(t2/tDelta); 
        nSampleTimeY = nSampleT1*2 + nSampleT2 + 1; caseID = caseID + 4; 
    end
    
    nSampleTime = max(nSampleTimeX,nSampleTimeY);
    switch caseID
        case 2
            t1 = tDelta*(nSampleTime-1)/2;
            if t1~=0
                ampGrad = flagKx*kx/(gamma/2/pi)/t1 + flagKy*ky/(gamma/2/pi)/t1;
                rateSlew = ampGrad/t1;
            end
        case 3
            nSampleTime = floor(nSampleTime/2)*2+1;
            t1y = tDelta*(nSampleTime-1)/2; 
            if t1y~=0
                ampGrady = ky/(gamma/2/pi)/t1y; rateSlewy = ampGrady/t1y;
            end                       
            nSampleT1x=nSampleTime-1-ceil(kx/(gamma/2/pi)/tDelta/ampGradMax);
            nSampleT2x = nSampleTime - 1 - nSampleT1x*2;
            t1x = nSampleT1x*tDelta; t2x = nSampleT2x*tDelta;
            ampGradx = kx/(gamma/2/pi*(t1x+t2x)); rateSlewx = ampGradx/t1x; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;              
        case 4
            nSampleTime = floor(nSampleTime/2)*2+1;
            t1x = tDelta*(nSampleTime-1)/2; 
            if t1x~=0
                ampGradx = kx/(gamma/2/pi)/t1x; rateSlewx = ampGradx/t1x;
            end
            nSampleT1y=nSampleTime-1-ceil(ky/(gamma/2/pi)/tDelta/ampGradMax);
            nSampleT2y = nSampleTime - 1 - nSampleT1y*2;
            t1y = nSampleT1y*tDelta; t2y = nSampleT2y*tDelta;
            ampGrady = ky/(gamma/2/pi*(t1y+t2y)); rateSlewy = ampGrady/t1y; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;     
        case 5
            nSampleT1x = ceil((nSampleTime-1)/2-sqrt((nSampleTime-1)^2-...
                4*kx/(gamma/2/pi)/rateSlewMax/(tDelta^2))/2);            
            nSampleT2x = nSampleTime - 1 - nSampleT1x*2;
            t1x = nSampleT1x*tDelta; t2x = nSampleT2x*tDelta;
            ampGradx = kx/(gamma/2/pi*(t1x+t2x)); rateSlewx = ampGradx/t1x; 
            nSampleT1y = ceil((nSampleTime-1)/2-sqrt((nSampleTime-1)^2-...
                4*ky/(gamma/2/pi)/rateSlewMax/(tDelta^2))/2);
            nSampleT2y = nSampleTime - 1 - nSampleT1y*2;
            t1y = nSampleT1y*tDelta; t2y = nSampleT2y*tDelta;
            ampGrady = ky/(gamma/2/pi*(t1y+t2y)); rateSlewy = ampGrady/t1y; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;   
        otherwise
    end
end

else

nSampleTimeRequired = round(tDuration/tDelta);
    
if imag(ka)*real(ka)==0
    flagKa = ka/abs(ka); ka = abs(ka);
    if ampGradMax^2/rateSlewMax >= ka/(gamma/2/pi)
        t1 = sqrt(ka/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1/tDelta); t1 = nSampleT1*tDelta;
        ampGrad = ka/(gamma/2/pi)/t1; rateSlew = ampGrad/t1; 
        nSampleTime = nSampleT1*2 + 1;
        if nSampleTimeRequired > nSampleTime
            if rem(nSampleTimeRequired,2) == 1
                nSampleT1 = floor((nSampleTimeRequired-1)/2);
                t1 = nSampleT1*tDelta;
                rateSlew = ka/(gamma/2/pi)/(t1^2 + t1*tDelta + tDelta^2);
                ampGrad = rateSlew * (t1+tDelta);
            else                
                nSampleT1 = nSampleTimeRequired/2; t1 = nSampleT1*tDelta;
                ampGrad = ka/(gamma/2/pi)/t1; rateSlew = ampGrad/t1; 
            end
        end
    else
        t1 = ampGradMax/rateSlewMax; t2 = ka/(gamma/2/pi)/ampGradMax-t1;
        nSampleT1 = ceil(t1/tDelta); t1 = nSampleT1*tDelta;
        nSampleT2 = ceil(t2/tDelta); t2 = nSampleT2*tDelta;
        ampGrad = ka/(gamma/2/pi*(t1+t2)); rateSlew = ampGrad/t1; 
        nSampleTime = nSampleT1*2 + nSampleT2 + 1;
        if nSampleTimeRequired > nSampleTime
            t1 = tDuration/2 - sqrt(tDuration.^2-4*ka/(gamma/2/pi)/rateSlewMax)/2;
            nSampleT1 = ceil(t1/tDelta); t1 = nSampleT1*tDelta;
            t2 = tDuration - t1*2; 
            nSampleT2 = ceil(t2/tDelta); t2 = nSampleT2*tDelta;            
            rateSlew = rateSlewMax; ampGrad = rateSlewMax*t1;
            kb = gamma/2/pi*ampGrad*(t1+t2);
        end
    end
    nSampleTime = nSampleTimeRequired;
    ampGrad = ampGrad*flagKa; rateSlew = rateSlew*flagKa; kb = kb*flagKa;
else    
    kx = real(ka); ky = imag(ka); caseID = 0;       
    flagKx = kx/abs(kx); kx = abs(kx); flagKy = 1i*ky/abs(ky); ky = abs(ky);
    if ampGradMax^2/rateSlewMax >= kx/(gamma/2/pi)
        t1x = sqrt(kx/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1x/tDelta); nSampleTimeX = nSampleT1*2 + 1;
        caseID = caseID + 0; 
    else
        t1x = ampGradMax/rateSlewMax; t2 = kx/(gamma/2/pi)/ampGradMax-t1x;
        nSampleT1 = ceil(t1x/tDelta); nSampleT2 = ceil(t2/tDelta); 
        nSampleTimeX = nSampleT1*2 + nSampleT2 + 1; caseID = caseID + 1; 
    end 
    if ampGradMax^2/rateSlewMax >= ky/(gamma/2/pi)
        t1y = sqrt(ky/(gamma/2/pi)/rateSlewMax);  
        nSampleT1 = ceil(t1y/tDelta); nSampleTimeY = nSampleT1*2 + 1;
        caseID = caseID + 2; 
    else
        t1y = ampGradMax/rateSlewMax; t2 = ky/(gamma/2/pi)/ampGradMax-t1y;
        nSampleT1 = ceil(t1y/tDelta); nSampleT2 = ceil(t2/tDelta); 
        nSampleTimeY = nSampleT1*2 + nSampleT2 + 1; caseID = caseID + 4; 
    end
    
    nSampleTime = max(nSampleTimeX,nSampleTimeY);
    switch caseID
        case 2
            t1 = tDelta*(nSampleTime-1)/2;
            if t1~=0
                ampGrad = flagKx*kx/(gamma/2/pi)/t1 + flagKy*ky/(gamma/2/pi)/t1;
                rateSlew = ampGrad/t1;
            end
        case 3
            nSampleTime = floor(nSampleTime/2)*2+1;
            t1y = tDelta*(nSampleTime-1)/2; 
            if t1y~=0
                ampGrady = ky/(gamma/2/pi)/t1y; rateSlewy = ampGrady/t1y;
            end                       
            nSampleT1x=nSampleTime-1-ceil(kx/(gamma/2/pi)/tDelta/ampGradMax);
            nSampleT2x = nSampleTime - 1 - nSampleT1x*2;
            t1x = nSampleT1x*tDelta; t2x = nSampleT2x*tDelta;
            ampGradx = kx/(gamma/2/pi*(t1x+t2x)); rateSlewx = ampGradx/t1x; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;              
        case 4
            nSampleTime = floor(nSampleTime/2)*2+1;
            t1x = tDelta*(nSampleTime-1)/2; 
            if t1x~=0
                ampGradx = kx/(gamma/2/pi)/t1x; rateSlewx = ampGradx/t1x;
            end
            nSampleT1y=nSampleTime-1-ceil(ky/(gamma/2/pi)/tDelta/ampGradMax);
            nSampleT2y = nSampleTime - 1 - nSampleT1y*2;
            t1y = nSampleT1y*tDelta; t2y = nSampleT2y*tDelta;
            ampGrady = ky/(gamma/2/pi*(t1y+t2y)); rateSlewy = ampGrady/t1y; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;     
        case 5
            nSampleT1x = ceil((nSampleTime-1)/2-sqrt((nSampleTime-1)^2-...
                4*kx/(gamma/2/pi)/rateSlewMax/(tDelta^2))/2);            
            nSampleT2x = nSampleTime - 1 - nSampleT1x*2;
            t1x = nSampleT1x*tDelta; t2x = nSampleT2x*tDelta;
            ampGradx = kx/(gamma/2/pi*(t1x+t2x)); rateSlewx = ampGradx/t1x; 
            nSampleT1y = ceil((nSampleTime-1)/2-sqrt((nSampleTime-1)^2-...
                4*ky/(gamma/2/pi)/rateSlewMax/(tDelta^2))/2);
            nSampleT2y = nSampleTime - 1 - nSampleT1y*2;
            t1y = nSampleT1y*tDelta; t2y = nSampleT2y*tDelta;
            ampGrady = ky/(gamma/2/pi*(t1y+t2y)); rateSlewy = ampGrady/t1y; 
            ampGrad = flagKx*ampGradx + flagKy*ampGrady; 
            rateSlew = flagKx*rateSlewx + flagKy*rateSlewy;   
        otherwise
    end
end    
    
end
% disp(ampGrad);disp(rateSlew);disp(nSampleTime);
end