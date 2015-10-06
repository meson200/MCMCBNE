function plotKSresult(CElist,CErand)

% plot the first KS filtering result (optimal dimensionality selection)

nrep = size(CElist,1);
nvar = size(CElist,2);
CE = mean(CElist,1);
CE_ste = 2*std(CElist,1)/sqrt(nrep-1);

CErand_1 = mean(CErand);
CErand_2 = 2*std(CErand)/sqrt(nrep-1);
rand_u = CErand_1+CErand_2;
rand_l = CErand_1-CErand_2;

figure
hold on
errorbar(1:length(CE),CE,CE_ste);
set(gca,'Xtick',1:length(CE))
set(gca,'FontSize',20)
xlabel('Number of variables eliminated','Fontsize',20);
ylabel('Cross entropy of the eliminated variable','Fontsize',20);
area([0 length(CE)],[rand_u rand_u],rand_l,'FaceColor','r','EdgeColor','None');
xlim([0 nvar+1]);
hold off