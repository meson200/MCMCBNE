function kyu_ReliabilityPlot(Perf_BN,Nbins,EnsSizesToPlot)

% reliability plot (estimated vs. actual probability)
% inputs:
% Perf_BN: prediction performance object from kyu_Perf_632BSplus
% Nbins: number of risk groups
% EnsSizesToPlot: an array of ensemble sizes to plot  (maximum 3)

Npatients = size(Perf_BN.Prob632p,1);
events = Perf_BN.Prob632p(:,2);
ens = Perf_BN.EnsembleSizes;
P_632p_BN = Perf_BN.Prob632p(:,3:end);
P_632p_BN_CI = Perf_BN.Prob632p_CI(:,3:end);
style = {'k*','ks','k^'};
Nmodels = numel(EnsSizesToPlot);
if Nmodels>3
    error('choose the number of sizes less than 3.')
end
sizes = [];
for i = 1:Nmodels
   sizes = [sizes find(ens==EnsSizesToPlot(i))];
end

legtext = repmat({'BN ensemble size'},[Nmodels 1]);
ensembles =  strread(num2str(EnsSizesToPlot),'%s');
mumatrix = P_632p_BN(:,sizes);
CImatrix = P_632p_BN_CI(:,sizes);

y = events-1;


figure
plot([0 1],[0 1],'rx--');
set(gca,'fontsize',20);
ylim([-0.01 1]);
xlim([0 1]);

[~,rx]=sort(mumatrix(:,end));


for p = 1:Nmodels
   
    mu = mumatrix(:,p);
    mu = abs(mu);
    CI_x = CImatrix(:,p);
    
    %[xs,rx]=sort(mu);
    %ys=y(rx);
    
    xs = mu(rx);
    ys = y(rx);

    logloss = 0;
    logloss_unc = 0;
    for i = 1:Npatients
        logloss = logloss - (ys(i)*log(xs(i)) + (1-ys(i))*log(1-xs(i)));
    end
    disp(['log loss = ',num2str(logloss)]);
    
    CI_xs = CI_x(rx);
    nBins=round(Npatients/Nbins);
    for i=1:Nbins
        
%         finish = start+vecsize(p,i)-1;
%         vec = start:finish;
%         start = finish+1;

        vec=(i-1)*nBins+[1:nBins];
        vec=vec(find(vec<=Npatients));
        
        yh(i)=mean(ys(vec));
        %if yh(i)==0
        %   yh(i) = yh(i) + 0.001; 
        %end
        %syh(i)=sqrt(yh(i)*(1-yh(i)))/sqrt(length(ys(vec)));
        syh(i) = sqrt(yh(i)*(1-yh(i))/length(vec));
        xh(i)=mean(xs(vec));
        sxh(i)=sqrt(xh(i)*(1-xh(i)))/sqrt(length(xs(vec)));
        rxh(i)=(yh(i)-xh(i))/sqrt(xh(i)*(1-xh(i)));
        sxh_x_tA = std(xs(vec))/sqrt(length(vec));
        sxh_x_tB = rms(CI_xs(vec));
        sxh_x(i) = sqrt(sxh_x_tA^2+sxh_x_tB^2);        
    end
    hold on
    h = errorbarxy(gca,xh,yh,sxh_x,syh,{style{p},'k','k'});
    yresid = yh-xh;
    SSresid = sum(yresid.^2);
    SStotal = (length(yh)-1) * var(yh);
    r2 = 1-SSresid/SStotal;
    r2str = strcat(' (',sprintf('%.2f',r2),')');
    set(h.hMain,'MarkerSize',15);
    h_ax(p) = h.hMain;
    lt = legtext{p};
    lt = [lt,' ',ensembles{p}];
    temp = strcat(lt,r2str);
    legtext{p} = temp;
end
legend(h_ax,legtext,'FontSize',20,'location','NorthWest');
xlabel('predicted risk','FontSize',20);
ylabel('actual risk','FontSize',20);
hold off



