function Wg = kyu_BN_FindRobustStruct(MCMCresults,data_in,ProbsObject)

% build a summarizing graph showing robust links and polarity of the links
% inputs
% this function can be called after a parameter training is finished
% MCMCresults: output of MCMC simulation (output of kyu_BN_GraphLearning_collect.m)
% data_in: folds of data used for graph/parameter training
% polarity: a matrix indicating polarity of correlation
% (output of kyu_BN_ParamLearnPredict_collect.m)
% polarity(i,j) = 1/-1: pos./neg. correlation between a variable i and j
% reference: Scutari (2013), On the Prior and Posterior Distributions Used
% in Graphical Modelling, Bayesian Analysis 8(1), pp.1-28

InputCorrect = false;
while ~InputCorrect
    resp = input('enter the ensemble size to average: ','s');
    try
        temp = textscan(resp,'%f');
        InputCorrect = ~isempty(temp);
        InputNo = temp{1};
    catch
        disp('input not numeric!');
        InputCorrect = false; 
    end
end
[~,ix] = min(abs(ProbsObject.EnsembleSizes-InputNo));
NoGraphs = ProbsObject.EnsembleSizes(ix);
disp(['The closest match: an ensemble size ',num2str(NoGraphs)]);

gs_top20 = MCMCresults.gs_top(:,1:NoGraphs);
post = MCMCresults.posterior(:,1:NoGraphs);
labels = data_in.KM.Labels;
polarity = ProbsObject.polarity{ix};

NumFolds = size(gs_top20,1);
NumStruct = size(gs_top20,2);
NumNodes = size(gs_top20{1},1);
Wg = zeros(size(gs_top20{1,1}));
Cnt = zeros(size(gs_top20{1,1}));
% determine a cutoff value on link confidence (see Scutari (2013))
% as a function of the number of nodes
conf_cutoff = 0.25+1/(4*(NumNodes-1)); 
% take an average of BN graphs in an ensemble to yield a confidence level at each link 
% when taking an average, each graph is weighted by posterior
WgCell = repmat({zeros(1,NumFolds)},NumNodes,NumNodes);
for i = 1:NumFolds
    Wgpart = 0;
    for j = 1:NumStruct
        newmat = gs_top20{i,j};
        Wgpart = Wgpart + newmat*post(i,j);
    end
    Wgpart = Wgpart/sum(post(i,1:NumStruct));
    for d1 = 1:NumNodes
        for d2 = 1:NumNodes
            WgCell{d1,d2}(i) = Wgpart(d1,d2);
        end
    end
end
% convert the average graph to a robust structure 
% using the cutoff value
Wg = cellfun(@mean,WgCell);
for i = 1:size(Wg,1)
    for j = 1:size(Wg,1)
        if Wg(i,j)>conf_cutoff && ttest(WgCell{i,j}-conf_cutoff)
            Cnt(i,j) = 1; 
        end
    end
end
% resolves the case where two variables are connected in both directions
% the link with lower confidence 
for i = 1:size(Wg,1)
    for j = 1:size(Wg,2)
        if Cnt(i,j) == 1 && Cnt(j,i) == 1
            if Wg(i,j)>Wg(j,i)
                Cnt(j,i) = 0;
            end
 
        end
    end
end
display which links were above the confidence cutoff
WgSparse = sparse(Wg);
[row,col] = find(WgSparse);
Ws = nonzeros(WgSparse);
[WsSorted,I] = sort(Ws,'descend');
for i = 1:length(I)
    indx = row(I(i));
    indy = col(I(i));
    if Ws(I(i))>conf_cutoff && ttest(WgCell{indx,indy}-conf_cutoff)
       str = sprintf('%s -> %s : delta=%1.4f, polarity=%d',labels{indx},labels{indy},Wg(indx,indy),sign(polarity(indx,indy)));
       disp(str);   
    end
end

% plot the robust structure using the MATLAB biograph 
connected = [];
for j = 1:size(Wg,1)
   if ~(1 && all(Cnt(j,:)==0) && all(Cnt(:,j)==0)) 
        connected = [connected j];
   end    
end
Cnt_s = Cnt(connected,connected);
Wg_s = Wg(connected,connected);
polarity_s = polarity(connected,connected);
labell = labels(connected);   

[par,chi] = find(Cnt_s);
bg = biograph(Cnt_s, labell,'Scale',2,'arrowsize',16);  

for i = 1:numel(find(Cnt_s))
    edgehandle = getedgesbynodeid(bg,labell{par(i)},labell{chi(i)});
    % thickness of an edge varies with a confidence level
    set(edgehandle, 'LineWidth', fn_edgethick(Wg_s(par(i),chi(i)),conf_cutoff));
    if polarity_s(par(i),chi(i)) > 0 % promotional relationship: blue
        ColorVec = [0 0 1];
        LinePattern = '-';
    else % inhibitory relationship: red
        ColorVec = [1 0 0];
        LinePattern = '-';
    end
    set(edgehandle, 'LineColor', ColorVec);
end
set(bg.Nodes,'FontSize',25);
%legend('P>0.95', 'Location','SouthWest', 'Color','g');
view(bg);

% function that maps a confidence level to an edge thickness displayed in
% the biograph
function width = fn_edgethick(weight,cutoff)
% if weight > 0.95 
%     width = 6;
% elseif weight > 0.68    
%     width = 3;
% elseif weight > conf_cutoff
%    width = 2;
% elseif weight > conf_cutoff
%     width = 1;
% end

maxW = 1;
minW = cutoff;
maxT = 9;
minT = 1;
width = (weight-minW)*(maxT-minT)/(maxW-minW);

% temporary function made to enable double indexing
function listout = SetValue(listin,index,value)

listout = listin;
listout(index) = value;





