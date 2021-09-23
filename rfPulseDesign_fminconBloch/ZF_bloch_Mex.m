function [rmse,mz] = ZF_bloch_Mex(pulseRf,mapB1,mapB0,spatial,targetVecDeg,gr,tDelta)
    %Convert the current rf pulse back from the polar form to the rectangular form
    pulseRf = pulseRf(1:end/2).*exp(1i*pulseRf(end/2+1:end));
    pulseRf = reshape(pulseRf,length(pulseRf)/8,8);
    mag = bloch3D_mex(mapB0,spatial,pulseRf,mapB1,gr,tDelta);
    mz = mag(3,:).';
    rmse = targetVecDeg - acosd(mz);
    rmse = sqrt(sum(rmse.^2)/length(targetVecDeg));
end