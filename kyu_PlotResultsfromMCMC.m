function kyu_PlotResultsfromMCMC(R,timept)

% plot posterior distribution from MCMC graph searching
% inputs:
% R: MCMC result object (created by kyu_BN_MCMC_convtest.m)
% timept: 6x1 array of iteration numbers to plot
% (should choose the numbers from the SampledT_large)


plotend = 500; % plot the highest posterior graphs 1~plotend
x = 1:plotend;

% identify which histograms to plot
t = zeros(6,1);
lstr = cell(6,1);
for i = 1:6
    t(i) = find(R.SampledT_large==timept(i));
    % for the legend
    lstr{i} = strcat('T=',num2str(timept(i)));
end

hist1 = R.DAGhistogram(t(1)).frequency;
hist2 = R.DAGhistogram(t(2)).frequency;
hist3 = R.DAGhistogram(t(3)).frequency;
hist4 = R.DAGhistogram(t(4)).frequency;
hist5 = R.DAGhistogram(t(5)).frequency;
hist6 = R.DAGhistogram(t(6)).frequency;

% generate cumulative histogram
histcum1 = zeros(1,500);
histcum2 = zeros(1,500);
histcum3 = zeros(1,500);
histcum4 = zeros(1,500);
histcum5 = zeros(1,500);
histcum6 = zeros(1,500);

for i = 1:plotend
    histcum1(i) = sum(hist1(1:i));
    histcum2(i) = sum(hist2(1:i));
    histcum3(i) = sum(hist3(1:i));
    histcum4(i) = sum(hist4(1:i));
    histcum5(i) = sum(hist5(1:i));
    histcum6(i) = sum(hist6(1:i));
end

% plot the posterior

hFig = figure(1);
set(hFig, 'Position', [1 1 1300 500])

subplot(1,2,1)
p1 = plot(x,hist6(1:plotend),'-','linewidth',2);
set(gca,'fontsize',20);
hold on
xlabel('Graph space','Fontsize',20)
ylabel('Posterior','Fontsize',20)
p2 = plot(x,hist2(1:plotend),':','linewidth',2);
p3 = plot(x,hist3(1:plotend),'--','linewidth',2);
p4 = plot(x,hist4(1:plotend),'b-',x(1:plotend/10:plotend),hist4(1:plotend/10:plotend),'x','linewidth',2);
p5 = plot(x,hist5(1:plotend),'b-',x(1:plotend/10:plotend),hist5(1:plotend/10:plotend),'o','linewidth',2);
p6 = plot(x,hist6(1:plotend),'b-',x(1:plotend/10:plotend),hist6(1:plotend/10:plotend),'d','linewidth',2);
leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],lstr{1},lstr{2},lstr{3},lstr{4},lstr{5},lstr{6});
%leg1 = legend([p1,p2,p3,p4(2),p5(2),p6(2)],'T=10000','T=20000','T=30000','T=40000','T=50000','T=60000');
set(leg1,'FontName','Helvetica');
pbaspect([2 2 1])
hold off

% plot the cumulative posterior
subplot(1,2,2)
axis square
p1 = plot(histcum1(1:plotend),'-','linewidth',2);
set(gca,'fontsize',20);
hold on
xlabel('Graph space','Fontsize',20)
ylabel('Posterior (cumulative)','Fontsize',20)
p2 = plot(x,histcum2(1:plotend),':','linewidth',2);
p3 = plot(x,histcum3(1:plotend),'--','linewidth',2);
p4 = plot(x,histcum4(1:plotend),'b-',x(1:plotend/10:plotend),histcum4(1:plotend/10:plotend),'x','linewidth',2);
p5 = plot(x,histcum5(1:plotend),'b-',x(1:plotend/10:plotend),histcum5(1:plotend/10:plotend),'o','linewidth',2);
p6 = plot(x,histcum6(1:plotend),'b-',x(1:plotend/10:plotend),histcum6(1:plotend/10:plotend),'d','linewidth',2);
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
