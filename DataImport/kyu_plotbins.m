function kyu_plotbins(data_C,data_C_raw,boundary,labels,units)

nvar = size(data_C,2);
plot_ncols = 8;
plot_nrows = ceil(size(data_C,2)/plot_ncols);
nbins = 7;
bins = linspace(-2,2,nbins);

for i = 1:nvar
    %subplot(plot_nrows,plot_ncols,i);
    subaxis(plot_nrows,plot_ncols,i,'sh',0.02,'sv',0.12,'Padding',0,'mt',0.08,'ml',0.03,'mr',0.03);
    isbinary = numel(unique(data_C_raw(:,i)))==2;
    if ~isbinary
        [nb,xb] = hist(data_C(:,i),bins);
        tt = data_C(:,i) - boundary(i);
        tt(tt<0) = 10^15;
        [~,b_raw_index] = min(tt);
    else
        temp = 2*(data_C_raw(:,i)-1.5);
        [nb,xb] = hist(temp,bins);
        boundary(i) = boundary(i) - (xb(2)-xb(1))/2;
    end
    wd = xb(2)-xb(1);
    hold on
    jlist = [];
    for j = 1:nbins
        bh=bar(xb(j),nb(j),wd);
        if xb(j)-wd/2 > boundary(i)
            set(bh,'facecolor',[1 0 0]);
            jlist = [jlist j];
        else
            set(bh,'facecolor',[0 0 1]);
        end
        
    end
    
    if strcmp(labels{i},'PTVCOMSI')    
        tickval = data_C_raw(b_raw_index,i);
        set(gca,'fontsize',20,'XTick',xb(min(jlist)),'XTicklabel',tickval);
        title([labels{i}],'Interpreter','none');
    elseif ~isbinary
        tickval = data_C_raw(b_raw_index,i);
        set(gca,'fontsize',20,'XTick',xb(min(jlist)),'XTicklabel',tickval);
        title([labels{i},'(',units{i},')'],'Interpreter','none');
    else
        set(gca,'fontsize',20,'XTick',-0.7,'XTicklabel',units{i});
        title([labels{i}],'Interpreter','none');
    end
    xlim([min(xb)-1 max(xb)+1])
    hold off

    
end


