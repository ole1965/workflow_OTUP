%% Evaluation of the results produced by 'main_databaseSizeInvestigation.m'.
%set the path where the results are stored.
path = 'databaseSizeInvestigation\results';
files = dir(path);
r = [];rvv = [];
for i = 4:size(files,1)
    load([path,files(i).name],'rmse_speicher_parfor','rv');
    r = [r;cell2mat(rmse_speicher_parfor)];
    rvv = [rvv,rv];
end
l = size(r,2);
%plot rmses and mean rmse per database size
figure;plot(r(:,1:l)',':');hold on;
plot(mean(r(:,1:l)),'color','blue','linewidth',5);hold off;

%Percentile Method for CIs
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
l = size(r,2)-1;
ms = [(r(:,1)-r(:,l)),(r(:,2)-r(:,l)),(r(:,3)-r(:,l)),(r(:,4)-r(:,l)),(r(:,5)-r(:,l)),(r(:,6)-r(:,l)),(r(:,7)-r(:,l)),...
      (r(:,8)-r(:,l))];

CI = [];
for i = 1:l
    CI(i,:) = CIFcn(ms(:,i),95);
end
figure;plot(ms',':');
hold on
plot(1:l,CI(:,1),'LineWidth',3,'Color','red')
plot(1:l,CI(:,2),'LineWidth',3,'Color','red');%sum(CI(8:end,2))
plot(mean(ms),'LineWidth',3,'Color','blue')
hold off;
