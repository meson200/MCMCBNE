function kyu_PlotResultsfromMCMC(R,timept)

% plot posterior distribution from MCMC graph searching
% inputs:
% R: MCMC result object (created by kyu_BN_MCMC_convtest.m)
% timept: 6x1 array of iteration numbers to plot
% (should choose the numbers from the SampledT_large)


plotend = 1000; % plot the highest posterior graphs 1~plotend
x = 1:plotend;
n_t = numel(timept);

% identify which histograms to plot
t = zeros(n_t,1);
lstr = cell(n_t,1);
for i = 1:n_t
    t(i) = find(R.SampledT_large==timept(i));
    % for the legend
    lstr{i} = strcat('T=',num2str(timept(i)));
end

hist = cell(n_t,1);
histcum = cell(n_t,1);

for i = 1:n_t
    hist{i} = zeros(plotend,1);
    histcum{i} = zeros(plotend,1);
    rawhist = R.DAGhistogram(t(i)).frequency;
    tofill = min(plotend,length(rawhist));
    hist{i}(1:tofill) = rawhist(1:tofill);
    for j = 1:plotend
        histcum{i}(j) = sum(hist{i}(1:j));
    end
end
% plot the posterior

hFig = figure(1);
set(hFig, 'Position', [1 1 1300 500])
%plotpattern = {'-',':','--','-','-',''};
%dotpattern = 
subplot(1,2,1)
set(gca,'fontsize',20);
p1 = plot(x,hist{1}(1:plotend),':','linewidth',2);
hold on
xlabel('Graph space','Fontsize',20)
ylabel('Posterior','Fontsize',20)
p2 = plot(x,hist{2}(1:plotend),':','linewidth',2);
p3 = plot(x,hist{3}(1:plotend),'--','linewidth',2);
p4 = plot(x,hist{4}(1:plotend),'b-',x(1:plotend/10:plotend),hist{4}(1:plotend/10:plotend),'x','linewidth',2);
p5 = plot(x,hist{5}(1:plotend),'b-',x(1:plotend/10:plotend),hist{5}(1:plotend/10:plotend),'o','linewidth',2);
p6 = plot(x,hist{6}(1:plotend),'b-',x(1:plotend/10:plotend),hist{6}(1:plotend/10:plotend),'d','linewidth',2);
leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],lstr{1},lstr{2},lstr{3},lstr{4},lstr{5},lstr{6});
%leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],'T=10000','T=20000','T=30000','T=40000','T=50000','T=60000');
set(leg1,'FontName','Helvetica');
pbaspect([2 2 1])
hold off

% plot the cumulative posterior
subplot(1,2,2)
axis square
p1 = plot(histcum{1}(1:plotend),'-','linewidth',2);
set(gca,'fontsize',20);
hold on
xlabel('Graph space','Fontsize',20)
ylabel('Posterior (cumulative)','Fontsize',20)
p2 = plot(x,histcum{2}(1:plotend),':','linewidth',2);
p3 = plot(x,histcum{3}(1:plotend),'--','linewidth',2);
p4 = plot(x,histcum{4}(1:plotend),'b-',x(1:plotend/10:plotend),histcum{4}(1:plotend/10:plotend),'x','linewidth',2);
p5 = plot(x,histcum{5}(1:plotend),'b-',x(1:plotend/10:plotend),histcum{5}(1:plotend/10:plotend),'o','linewidth',2);
p6 = plot(x,histcum{6}(1:plotend),'b-',x(1:plotend/10:plotend),histcum{6}(1:plotend/10:plotend),'d','linewidth',2);
%leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],'T=10000','T=20000','T=30000','T=40000','T=50000','T=60000');
leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],lstr{1},lstr{2},lstr{3},lstr{4},lstr{5},lstr{6});
set(leg1,'FontName','Helvetica','location','northwest');
pbaspect([2 2 1])
hold off

% Plot likelihood score & acceptance rate ( and beta)
%MCMC_length = length(R.score); 
PlotLength = max(timept);
ScoreAxisMax = max(R.score)+(max(R.score)-min(R.score))/2;
ScoreAxisMin = min(R.score)-(max(R.score)-min(R.score))/2;
ACRAxisMax = 0.8; % plotting range for acceptance rate

MaxScoreIndex = find(PlotLength == R.SampledT_small);
ACRm = zeros(1,MaxScoreIndex);
for i = 1:MaxScoreIndex
    ACRm(i) = sum(R.ACR((i-1)*10+1:i*10));
end
ACRm = ACRm/10;
Xaxis = R.SampledT_small(1:MaxScoreIndex);

figure
[ax h1 h2]=plotyy([0,Xaxis(end)],[ScoreAxisMin,ScoreAxisMax],[0,Xaxis(end)],[0,ACRAxisMax]); 
set(ax,'XScale','log');
cla(ax(1))
cla(ax(2))
axes(ax(1));
hold(ax(1))
p1 = plot(Xaxis,R.score,'-r','linewidth',2);
set(ax(1),'xticklab',[],'xtick',[],'FontSize',20,'YLim',[ScoreAxisMin,ScoreAxisMax],'YColor','k');
ylabel('Likelihood');
%set(ax(1),'YLim',[-300 -240],'Fontsize',20,'YLabel','Likelihood','YColor','-k');
%ylabel('Posterior','Fontsize',20)
%ylim([-300 -240]);
axes(ax(2));
hold(ax(2))
p2 = plot(Xaxis,ACRm,'-g','linewidth',2);
set(ax(2),'FontSize',20,'YColor','k','YLim',[0 ACRAxisMax])
ylabel('Acceptance rate');
%set(gca,'FontSize',20,'YColor','k','YLabel',20)
leg1 = legend([p1,p2],'score','acceptance rate','location','northwest');
xlabel('#runs','Fontsize',20)
hold off
