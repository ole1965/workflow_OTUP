%Script to perform the database size investigation as described in the manuscript and step 2 in the 'readMe.md'.
%in each iteration of the for loop a new random order of heads will be created.
for u = 1:500
    setUp;
    
    %the following four lines only work on a linux system
    delete myrand1.txt
    command = ['shuf -i 1-',num2str(nrDatasets ),' -n ',num2str(nrDatasets),' >> myrand1.txt'];
    system(command);
    rv = importdata('myrand1.txt');

    
    load(whichGrad)
    %pulse design for each database size
    parfor z = 1:nrDatasets 
    
        headsDB = 1:z;
        headsDB = rv(headsDB);
        headsTest = rv(end);

        %create design database
        mapB1 = []; mapB0 = []; mask = []; spatial = [];
        for j = 1:length(headsDB)
            switch target
                case 1 %wB
                    [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_wB(['datasets\recon_subj',num2str(headsDB(j)),'.mat']);
                case 2 %LEx
                    [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['datasets\recon_subj',num2str(headsDB(j)),'.mat']);
                    %postprocess the design mask according to the nuclei application in the paper
                    masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
                    masksingle(~mask2) = false;
            end
            masksingle = logical(masksingle(:));
            mapB1 = [mapB1;mapB1single];
            mapB0 = [mapB0;mapB0single];
            mask = [mask;masksingle];
            spatial = [spatial;spatialsingle];
        end

        maxiter = 20;
        [pulseRf,err,sta] = pulseRFcalcMLsqr(mapB1,mapB0,spatial,targetVec,gr,tDelta,maxiter);

        Gr = [real(gr(:,1)),imag(gr(:,1)),gr(:,2)]';
        %test the resulted pulse on the non-database head
        for j = 1:length(headsTest)
            disp(headsTest(j))
            switch target
                case 1 %wB
                    [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_wB(['datasets\recon_subj',num2str(headsTest(j)),'.mat']);
                case 2 %LEx
                    [mapB1single,mapB0single,masksingle,spatialsingle] = prepare_b1_b0_mask_LEx(['datasets\recon_subj',num2str(headsTest(j)),'.mat']);
                    %postprocess the design mask according to the nuclei application in the paper
                    masksingle = reshape(masksingle,[64 64 64]);masksingle(:,:,[1:startPos-1,endPos+1:end]) = 0;
                    masksingle(~mask2) = false;
            end
            masksingle = logical(masksingle(:));
            
            %Prepare the Data for the Mex Bloch Simulation
            mapB0single = mapB0single(logical(masksingle))';
            spatial = [real(spatialsingle(logical(masksingle),1)),imag(spatialsingle(logical(masksingle),1)),spatialsingle(logical(masksingle),2)]'; %dist = dist(:,logical(mskTarget));

            mapB1single = mapB1single(logical(masksingle),:);

            mag = bloch3D_mex(mapB0single,spatial,pulseRf,mapB1single,Gr,tDelta);
            mz = mag(3,:).';
            mzDeg = acosd(mz);

            targetVecSingleMasked = targetVecSingle; targetVecSingleMasked(~masksingle) = [];
            rmse = rad2deg(targetVecSingleMasked) - mzDeg;
            rmse = sqrt(sum(rmse.^2)/sum(masksingle));
            rmse_save{z} = rmse;
        end
    end
    %save the resulted rmse values
    t = datetime('now','Format','yyyy-MM-dd''T''HHmmssSSS');S = char(t);
    save(['databaseSizeInvestigation\results\results_',name,'_',S,'.mat'],'rmse_save','rv','whichGrad');
end