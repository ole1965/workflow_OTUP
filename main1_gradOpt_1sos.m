%Script to optimize the 1sos basis trajectory using the particle swarm optimization (as described in the manuscript and step 1 in the 'readMe.md').
setUp;
addpath('gradOpt');

switch target
    case 1 %targetWB
        [mapB1,mapB0,mask,spatial] = prepare_b1_b0_mask_wB(['B0B1datasets\recon_subj1.mat']);
        if startPos > 1 && endPos < 64
                mask = reshape(mask,[64 64 64]);mask(:,:,[1:startPos-1,endPos+1:end]) = 0;
        end
    case 2 %targetNuclei
        [mapB1,mapB0,mask,spatial] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj1.mat']);
        %postprocess the design mask according to the nuclei application in the paper
        if startPos > 1 && endPos < 64
            mask = reshape(mask,[64 64 64]);mask(:,:,[1:startPos-1,endPos+1:end]) = 0;
            mask(~mask2) = false;
        end        
    case 3 %targetM
        [mapB1,mapB0,mask,spatial] = prepare_b1_b0_mask_LEx(['B0B1datasets\recon_subj1.mat']);
        if startPos > 1 && endPos < 64
            mask = reshape(mask,[64 64 64]);mask(:,:,[1:startPos-1,endPos+1:end]) = 0;
        end
end
mask = logical(mask(:));

mapB1 = mapB1(logical(mask),:);
mapB0 = mapB0(logical(mask),:);
spatial = spatial(logical(mask),:);
targetVec = targetVec(logical(mask),:);

% parpool([5,128]);

%lb and ub depicts the lower and upper bounds for the spiral parameters in Eq. 1 in the manuscript
lb = [0.001 , 0, 0.001 , 3   , -100]; %[gamma,alpha,n (number of spiral turns, integrated in omega),spiral length,location of the spiral on the z axis]
ub = [200   , 5,    15 , 950 ,  100];

options = optimoptions('particleswarm','UseParallel',true,'Display','iter');

maxiter = 20; %number of iterations for the lsqr-algorithm inside the cost function ZF_1sos
nr = 20;
rmse_val = zeros(1,nr);paras = zeros(nr,length(lb));
%call the particleswarm optimization nr times
%paras is the vector to be optimized
for k = 1:nr
    [paras(k,:),rmse_val(k)] = particleswarm(@(paras) ZF_1sos(paras,mapB1,spatial,targetVec,maxiter,mapB0),...
                5,lb,ub,options);
    save gradientWaveforms\optParas\results_gradOpt_1sos.mat paras rmse_val
end

rmse_val(rmse_val == 0) = 1000;
[~,pos] = min(rmse_val); %Always check if there is another set of parameters that is also suitable,
%for instance because it provides a shorter pulse duration with simular rmse performance
paras = paras(pos,:);
gr = transformParasToWaveform_1sos(paras);

t = datetime('now','Format','yyyy-MM-dd''T''HHmmssSSS');S = char(t);
save(['gradientWaveforms\gradient_1sos_',targetstr,'_',S,'.mat'],'gr');