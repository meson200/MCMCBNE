function [mat_freq_out,hist_peaked,mat_top] = matrix_hist(mats,freq)

% generate posterior probability distribution from graph samples
% the posterior will be evaluated every T samples = time interval to
% assess convergence
%
% inputs
% mats: MCMC sampled graphs 
% freq: time interval to estimate posterior
% 
% outputs
% mat_freq_out: an array of a struct variable that has 2 components
% .matrix: matrix represenation of the graph
% .frequency: occurrence of the graph in the samples (max:1)
% length(mat_freq_out) = number of different graphs found in the MCMC samples
% hist_peaked: peakedness of the posterior graph. 
% defined as Ngraphs(0.5*highest frequency)/Ngraphs(mat_freq_out) 
% mat_top: a matrix 

PosteriorHist = {};
mat_freq_out = [];
mat_top = [];
hist_peaked = [];
n = length(mats);
mat_freq = struct('matrix',{},'frequency',{});
found = false;

p = 0;
for i = 1:n
    sz = size(mat_freq,2);    
    newmat = mats{i};
    if sz > 0
        found = false;
        for k = 1:sz
            if isequal(mat_freq(k).matrix,newmat) % if the matrix has been already found
                mat_freq(k).frequency = mat_freq(k).frequency + 1;
                found = true;
            end
        end
    end

    if sz == 0 || ~found % new dag entry
        mat_freq(sz+1).matrix = newmat;
        mat_freq(sz+1).frequency = 1;
    end
    
    if mod(i,freq) == 0
       p = p+1;
       histo = [];
       for m = 1:length(mat_freq)
            histo = [histo mat_freq(m).frequency];
       end
       [hist_sorted,sortedindex] = sort(histo,'descend');
       hist_sorted = hist_sorted/sum(hist_sorted);
       mat_freq_out(p).frequency = hist_sorted;
       PosteriorHist{p} = hist_sorted;
       for l = 1:numel(sortedindex)
           mat_freq_out(p).matrix{l} = mat_freq(sortedindex(l)).matrix; 
       end
       mat_top{p} = mat_freq_out(p).matrix{1};
       freqmax = hist_sorted(1);
       hist_peaked(p) = max(find(hist_sorted-freqmax/2>0))/numel(histo);
    end

end


%max2index = sortedindex(find(hist==));


% draw the top 
%bh = kyu_biograph(mat_top,labels);
%set(bh.Nodes,'FontSize',20);
%view(bh)
% optional: draw the second & third
% mat_second = mat_freq(maxindex2).matrix;
% bh2 = kyu_biograph(mat_second,labels);
% set(bh2.Nodes,'FontSize',20);
% view(bh2)
% mat_third = mat_freq(maxindex3).matrix;
% bh3 = kyu_biograph(mat_third,labels);
% set(bh3.Nodes,'FontSize',20);
% view(bh3)




% normalize the histogram
hist_sorted = hist_sorted/sum(hist_sorted);
