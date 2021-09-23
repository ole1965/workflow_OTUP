function [targetMat,mask2] = prepare_target(target,startPos,endPos,nrSl)
    switch target
        case 1 %targetWB
            targetMat = ones(64,64,64);
            mask2 = [];
        case 2 %targetNuclei
            load('targets\target2D_LEx_nuclei.mat');
            %prepare second mask to exclude the voxels directly adjacent to the excited area
            mask2 = ones(64,64);
            mask2(27:41,25:39) = 0;
            mask2 = mask2+target2D;
            mask2(mask2~=0) = 1;
            mask2 = logical(repmat(mask2,[1,1,nrSl]));
            mask2 = cat(3,zeros(64,64,startPos-1),mask2,zeros(64,64,64-endPos));
            
            targetMat = repmat(target2D,[1,1,nrSl]);
            targetMat = cat(3,zeros(64,64,startPos-1),targetMat,zeros(64,64,64-endPos));
        case 3 %targetM
            load('targets\target2D_LEx_M.mat');
            targetMat = repmat(target2D,[1,1,nrSl]);
            targetMat = cat(3,zeros(64,64,startPos-1),targetMat,zeros(64,64,64-endPos));
            mask2 = [];
    end
end