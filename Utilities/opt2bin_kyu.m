function [mi boundary binneddata] = opt2bin_kyu (rawdata, class, steps,...
    typesearch, minint, maxint)
% Karl Kuschner, College of William and Mary, Dept. of Physics, 2009.
% edited by Sangkyu Lee 2015
%
% opt2bin finds the best single boundary for each variable to maximize MI 
% 
% DESCRIPTION
%       This function takes an array of continuous data, with cases in rows
%       and variables in columns, along with a vector "class" which holds
%       the known class of each of the cases, and returns an array
%       "binneddata" that holds the 2 bin discretized data.  The
%       discretization bin boundary is found by maximizing the mutual
%       information with the class; the resulting MI and boundary are also
%       returned. The starting boundaries for the search can be given in
%       the vectors min and max, or either one, or neither, in which case
%       the data values determine the search boundaries.% 
%
% USAGE
%       [mi boundary binneddata] = maxMIbin(rawdata, class, typesearch [,
%           min, max])
% 
% INPUTS
%       rawdata: double array of continuous values, cases in rows and 
%           variables in columns. Distribution is unknown.
%       class: double column vector, values 1:c representing classification
%           of each case. 
%       steps: Number of steps to test at while finding maximum MI
%       typesearch =0: starting bndry based on data's actual max/min values
%                  =1: use the value passed in max as maximum (right) value
%                  =-1: use the value passed in min as minimum (left) value
%                  =2: used values passed via max, min
%       the two optional arguments are vectors whose values limit the range
%       of search for each variables boundaries.
% 
% OUTPUTS
%
%       mi: row vector holding the maximum values of MI(C;Vi) found
%       boundary: The location used to bin the data to get max MI
%       binneddata: The resulting data binned into "1" (low) or "2" (hi)
% 
% CALLED FUNCTIONS
% 
%       MIarray: Finds the MI of each col in an array with a separate
%           vector (the class in this case)

%% Intialize
[rows cols]=size(rawdata);
mi=zeros(1,cols);
boundary=zeros(1,cols);
binneddata=zeros(rows,cols);
currentmi=zeros(steps,cols);

% if not passed, find the left and rightmost possible bin boundaries from
% data

if nargin~=6
    minint=min(rawdata,[],1);
    maxint=max(rawdata,[],1);
elseif typesearch==1
    minint=min(rawdata,[],1);
elseif typesearch==-1
    maxint=max(rawdata,[],1);
elseif typesearch==2
    disp('using passed values')
else
    disp('typesearch must = 0,1,-1,2')
    return
end

%% Find best boundary

for peak=1:cols %look at each variable separately

    % Create an array of bin boundary's possible locations min->max
    checkpoints=repmat(linspace(minint(peak),maxint(peak),steps),rows,1); 
    
    % discretize the variable's values at each of these possible
    % boundaries, putting 2's everywhere (value > boundary), 1 elsewhere 
    binarray=(repmat(rawdata(:,peak), 1, steps)>checkpoints)+1; 
    
    % Send this array off to find the MI(C,V) for each possible binning  
    currentmi(1:steps,peak)=MIarray_kyu(binarray,class);


    % Now pick out the highest MI, i.e. best bin boundary 
    [mi(peak) atstep]=max(currentmi(:,peak));
    boundary(peak)=checkpoints(1,atstep);
    
    % and record the binned data using that boundary.
    binneddata(:,peak)=binarray(:,atstep);
end

end