function MIOut = MIarray(MatrixIn, class)
% by Karl Kuschner, College of William and Mary, Dept. of Physics, 2009.
%
% MIarray finds MI of each column of a data set with a separate vector
%
% DESCRIPTION
%       This function finds the mutual information between a single
%       discrete variable (class) and a matrix of discrete variables
%       (MatrixIn) which have the same number of cases (variables in
%       columns, cases in rows). A row vector containing the values
%       MI(Vi,C) for each variable Vi in the matrix is returned
%
% USAGE
%       MIOut = MI(MatrixIn, class)
%
% INPUTS
%       data: double array of discrete integer (1:n) values, cases in rows
%           and variables in columns.
%       class: double (col) vector, values 1:c representing class of each
%           case. Number of values c can be different than n in the data.
%
% OUTPUTS
%       MIOut: double (row) vector whose entries are the Mutual information
%           between each corresponding column of MatrixIn and the class.
%
% CALLED FUNCTIONS
%       None.
%


%% Intialize and Data Check
% check arguments
if nargin~=2
    disp('wrong number of input arguments')
    disp('need (data_array,class)')
    disp(' ')
    disp('Type "doc MI" for more info')
end
%class and MatrixIn must have the same number of rows
[rows cols]=size(MatrixIn);
if size(class,2)==rows
    class=transpose(class);
elseif size(class,1)~=rows
    disp('Dimension mismatch in rows of MI arrays')
    disp('Input arrays must have the same number of rows')
    return
end %row dimension check


% States must be integer values, typically 1 to n. If so, record n.
% Similarly, find out the number of states of the class variable.
if sum(any(MatrixIn-round(MatrixIn)))
    disp('Matrix in should be integers 1 to n')
    return
else
    n=max(size(unique(MatrixIn))); % Number of data states
    c=max(size(unique(class)));% Number of class states
end % check if integer

%% Variable Prep

MatrixIn=int8(MatrixIn); %optional
class=int8(class); %optional
Pcv = zeros(c,n,cols);

%% Main function

% Compute probability tables. P_ij is a matrix whose entries are
% Prob(Variable 1=state i and Variable 2= state j).  Others are similar.

if c==1 %trap for errors in the case where all classes are the same
    Pc = 1;
else
    % Create a 3-D array with c rows, each row filled with P(C=ci)
    Pc = repmat((hist(class,1:c)/rows)', [1,n,cols]);
end

% Create a 2-D array where (j,k) is P(Vk=vj).  Replicate it to a third
% dimension to prepare for multiplication with the above.
Pv =repmat( reshape ( hist(MatrixIn,1:n)/rows , [1,n,cols] ), [c,1,1] );

% Now multiply these together,  The result is a c by n by cols matrix whose
% (i,j,k) entry is P(C=ci)*P(Vk=vj) for each value of class ci and data vj.
PcPv= Pc.*Pv;

% Now we need a similar sized array with the  (i,j,k) entry equal to
% P(C=ci and Vk=vj) -- the joint probability.
for classstate=1:c
    Pcv(classstate,:,:) = hist(MatrixIn(class==classstate,:),1:n)/rows;
end

% Now we can compute the mutual info using
% 
% MI(C=i;Vk=j) = sum i (sum j (Pcv(i,j,k) log [Pcv(i,j,k)/PcPv(i,j,k)] ) )
% 
miterms=Pcv.*(log2(Pcv)-log2(PcPv)); % The term inside the log above...
miterms(isnan(miterms))=0; % with all the 0 log 0 entries removed

% Do the double summation and squeeze the unused dimensions
MIOut = squeeze(sum(sum(miterms,1),2))';

end