function imputed = knnimpute_kyu(data,K,varargin)
%KNNIMPUTE imputes missing data using the nearest-neighbor method
% edited from knnimpute.m by Kyu on May 5, 2015
% data is z-scored before calculating pairwise distance
%
%   KNNIMPUTE(DATA) replaces NaNs in DATA with the corresponding value from
%   the nearest-neighbor column using Euclidean distance. If the nearest-
%   neighbor column also contains a NaN value, then the next nearest column
%   is used.
%
%   KNNIMPUTE(DATA, K) replaces NaNs in DATA with a weighted mean of the K
%   nearest-neighbor columns. The weights are inversely proportional to the
%   distances from the neighboring columns.
%
%   KNNIMPUTE(..., 'DISTANCE',DISTFUN) computes nearest neighbors using
%   distance metric DISTFUN.  Choices are
%
%       'euclidean'   Euclidean distance
%       'seuclidean'  Standardized Euclidean distance, in which each
%                     coordinate in the sum of squares is inverse weighted
%                     by the sample variance of that coordinate
%       'cityblock'   City block distance
%       'mahalanobis' Mahalanobis distance
%       'minkowski'   Minkowski distance with exponent 2
%       'cosine'      One minus the cosine of the included angle
%                     between observations (treated as vectors)
%       'correlation' One minus the sample correlation between
%                     observations (treated as sequences of values)
%       'hamming'     Hamming distance, percentage of coordinates
%                     that differ
%       'jaccard'     One minus the Jaccard coefficient, the
%                     percentage of nonzero coordinates that differ
%       'chebychev'   Chebychev distance (maximum coordinate difference)
%       function      A distance function specified using @, for
%                     example @DISTFUN
%
%   See PDIST for more details.
%
%   KNNIMPUTE(..., 'DISTARGS',ARGS) passes the arguments ARGS to the
%   function DISTFUN. ARGS can be a single value or a cell array of
%   values.
%
%   KNNIMPUTE(...,'WEIGHTS',W) allows you to specify the weights used in
%   the weighted mean calculation. W should be a vector of length(K).
%
%   KNNIMPUTE(...,'MEDIAN',true) uses the median of the K nearest neighbors
%   instead of the weighted mean.
%
%   Examples:
%
%       load yeastdata
%       % Remove data for empty spots
%       emptySpots = strcmp('EMPTY',genes);
%       yeastvalues(emptySpots,:) = [];
%       genes(emptySpots) = [];
%       % Impute missing values
%       imputedValues = knnimpute(yeastvalues);
%
%   See also ISNAN, KNNCLASSIFY, NANMEAN, NANMEDIAN, PDIST.

% References: 
%   Speed, T. Statistical Analysis of Gene Expression Microarray Data
%   (2003), Chapman & Hall
%
%   Hastie, T., Tibshirani, R., Sherlock, G., Eisen, M., Brown, P. and
%   Botstein, D., Imputing Missing Data for Gene Expression Arrays, 
%   Stanford University Statistics Department Technical report (1999),
%   http://www-stat.stanford.edu/~hastie/Papers/missing.pdf
%
%   Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T.,
%   Tibshirani, R., Botstein, D.,  Altman, R.B., Missing value estimation
%   methods for DNA microarrays BIOINFORMATICS Vol. 17 no. 6 (2001) Pages
%   520-525 


% Copyright 2003-2009 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

if nargin < 2
    K = 1;
end
K = max(min(round(K),size(data,2)-1),1);

metric = 'euclidean';
uWeights = [];
userWeights = false;
distargs = {};
useWMean = true;

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:knnimpute:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'distance','distargs','median','weights'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        karg = find(strncmpi(pname, okargs,length(pname)));
        if isempty(karg)
            error(message('bioinfo:knnimpute:UnknownParameterName', pname));
        elseif length(karg)>1
            error(message('bioinfo:knnimpute:AmbiguousParameterName', pname));
        else
            switch(karg)
                case 1  % distance
                    metric = pval;

                case 2 %distfun args
                    if ~iscell(pval)
                        distargs = {pval};
                    else
                        distargs = pval;
                    end

                case 3 % median flag
                    medianFlag = bioinfoprivate.opttf(pval,okargs{karg},mfilename);
                    if medianFlag
                        useWMean = false;
                    end
                case 4  % weights
                    uWeights = pval(:);
                    userWeights = true;
            end
        end
    end
end

% create a copy of data for output
imputed = data;

% identify missing vals
nanVals = isnan(data);
% use rows without nans for calculation of nearest neighbors
noNans = sum(nanVals,2) == 0;

dataNoNans = data(noNans,:);
if isempty(dataNoNans)
    error(message('bioinfo:knnimpute:AllNans'));
end

% check weights
if userWeights
    if numel(uWeights) ~= K || ~isnumeric(uWeights)
        error(message('bioinfo:knnimpute:BadWeights'));
    end
    if ~useWMean
        warning(message('bioinfo:knnimpute:MedianWeights'));
    end
end

% For the special case of dependent columns in the matrix with non-nans, a
% more appropriate warning is given
if strncmpi(metric,'mahalanobis',max(length(metric),2)) && rank(cov(dataNoNans'))<size(dataNoNans,1)
    warning(message('bioinfo:knnimpute:DependentColumns'))
end

% normalize
dataNoNansNorm = zscore(dataNoNans,1,2);

% calculate pairwise distances between columns
distances = pdist(dataNoNans',metric,distargs{:});

% sort and get indices of nearest columns
SqF = squareform(distances);
% force the diagonals to be negative so that they automatically sort to the
% top
[dists, ndx] = sort(SqF - eye(size(SqF)));
% find locations where ordered neighbors are equal distanced
equalDists = [diff(dists(2:end,:),1,1)==0;false(1,size(dists,2))];

% get rows and columns of missing values
[rows,cols] = find(nanVals);

% for each missing value find the nearest column without a corresponding
% NaN.

rowWarn = false(size(data,1),1);
for count = 1:numel(rows)
    % check that we don't have a row of nans
    if all(isnan(data(rows(count),:)))
        if rowWarn(rows(count)) == false
            warning(message('bioinfo:knnimpute:RowAllNans', rows( count )));
            rowWarn(rows(count)) = true;
        end
        continue
    end
    % start at 2 as we know that each col is the closest element
    for nearest = 2:size(ndx,1)-K+1
        % Look for exactly equal columns to the last of the K neighbors and
        % add them into the calculation, L contains the number of extra
        % neighbors that need to be considered.
        L = find(equalDists(nearest+K-2:end,cols(count))==0,1)-1;
        
        dataVals = data(rows(count),ndx(nearest:nearest+K+L-1,cols(count)));
        if useWMean
            if ~userWeights
                weights = 1./dists(2:K+L+1,cols(count));
            else
                if L > 0
                    weights = [uWeights(1:end-1); repmat(uWeights(end)/(L+1),L+1,1)];
                else
                    weights = uWeights;
                end
            end
            val = wnanmean(dataVals,weights);
        else
            val = nanmedian(dataVals);
        end
        if ~isnan(val)
            imputed(rows(count),cols(count)) = val;
            break
        end
    end
end


function m = wnanmean(x,weights)
%WNANMEAN Weighted Mean value, ignoring NaNs, infs are special

% Find NaNs and set them to zero
x = x(:); weights = weights(:);
nans = isnan(x);
infs = isinf(weights);
if all(nans)
    m = NaN;
    return
end
if any(infs)
    m = nanmean(x(infs));
    return 
end
% Sum up non-NaNs, and divide by the number of non-NaNs.
x(nans) = 0;
weights(nans) = 0;
% normalize the weights
weights = weights./sum(weights);
m = weights'*x;
