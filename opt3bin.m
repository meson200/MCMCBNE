function [l, r, binned, mi] = opt3bin (data, class)
% by Karl Kuschner, College of William and Mary, Dept. of Physics, 2009.
%
% FunctionName short description 
% 
% DESCRIPTION
%       This function takes an array of continuous sample data of size
%       cases (rows) by variables (columns), along with a class vector of
%       integers 1:c, each integer specifying the class. The class vector 
%       has the same number of cases as the data.  The function outputs the
%       position of the 2 bin boundaries (3 bins) that optimize the mutual
%       information of each variable's data vector with the class vector.    
% 
% USAGE
%       [l,r,binned, mi]=opt3bin(data,class)
% 
% INPUTS
%       data: double array of continuous values, cases in rows and 
%           variables in columns. Distribution is unknown.
%       class: double column vector, values 1:c representing classification
%           of each case. 
% 
% OUTPUTS
%
%       l     - row vector of left boundary position for each var.
%       r     - row vector of right boundary position for each var.
%       binned- data array discretized using boundaries in l and r
%       mi    - row vector of mutual info between each discr. variable 
%                  and class 
% 
% CALLED FUNCTIONS
% 
%       opt2bin: Similar function that finds a single boundary. This is
%           used as a seed for the 3 bin optimization.
%       looklr: See below.


%% Intialize
% 
%  Variable Prep : find sizes of arrays and create placeholders for locals

%steps=150;
steps = 300;
[rows cols]=size(data);
boundary=zeros(2,cols);

%% Method
% Find starting point by finding the maximum value of a 2 bin mi. Next, go
% left and right from that position, finding the position of the
% next boundary that maximizes MI.

[mi boundary(1,:)] = opt2bin (data, class, steps, 2);

% We've located a good starting (center) bin boundary.  Search L/R for a
% second boundary to do a 3 bin discretization.
[mi boundary(2,:)] = looklr (data, class, boundary(1,:), steps);

% We've now found the optimum SECOND boundary position given the best 2 bin
% center boundary.  Now re-search using that SECOND boaundary position,
% dropping the original (2 bin).  The result should be at, or near, the
% optimal 3 bin position.
[mi boundary(1,:) binned] = looklr (data, class, boundary(2,:), steps);

% from the two boundaries found above, sort the left and right
r=max(boundary);
l=min(boundary);

% Now retutn the vector of left and right boundaries, the disc. data, and
% max MI found.
end % of function


function [miout nextboundary binned] = looklr (data, class, startbd, steps)
% given a start position, finds another boundary (to create 3 bins) that
% maximizes MI with the class
[rows cols]=size(data);
farleft=min(data,[],1);
farright=max(data,[],1);
miout=zeros(1,cols);
binned=zeros(rows,cols);
nextboundary=zeros(1,cols);

for peak=1:cols % for each peak/variable separately...

    % discretize this variables' values. Sweep through the possible
    % bin boundaries from the startbd to the furthest value of the
    % data, creating 2 boundaries for 3 bins. Record the binned values in
    % a "cases x steps" array, where "steps" is the granularity of the
    % sweep. The data vector starts off as a column...

    testmat=repmat(data(:,peak),1,steps); % and is replicated to an array.
    
    % Create same size array of bin boundaries. Each row is the same.
    checkptsL=repmat(linspace(farleft(peak),startbd(peak),steps),rows,1); 
    checkptsR=repmat(linspace(startbd(peak),farright(peak),steps),rows,1);
    
    % Create a place to hold the discrete info, starting with all ones. The
    % "left" array will represent data binned holding the center boundary
    % fixed and sweeping out a second boundary to the left; similarly the
    % right boundary starts at "startbd" and sweeps higher.
    binarrayL=ones(rows,steps); 
    binarrayR=ones(rows,steps); 
    
    % Those in the L test array that are higher than the left boundary -> 2
    binarrayL(testmat>checkptsL)=2;
    binarrayL(testmat>startbd(peak))=3; % >center boundary -> 3
    
    % Similarly using center and right boundaries
    binarrayR(testmat>startbd(peak))=2;
    binarrayR(testmat>checkptsR)=3;

    % Now at each of those step positions, check MI (var;class).
    miout(peak) = 0;
    
    % THese vectors hold the MI with each step used to discretize.
    miL=MIarray(binarrayL,class);% MI(V;C) using left/center
    miR=MIarray(binarrayR,class); % MI(V;C) using center/right

 if max(miL)>max(miR)  % See which one is the largest
     [miout(peak) index]=max(miL); %record the max mi found
     nextboundary(peak)=checkptsL(1,index); % and record the boundary
     binned(:,peak)=binarrayL(:,index);% and record the discrete data
 else
     [miout(peak) index]=max(miR); %record the max mi found
     nextboundary(peak)=checkptsR(1,index); % and record the boundary
     binned(:,peak)=binarrayR(:,index);% and record the discrete data
 end
    
end % of that variable's search.  Go to next variable.

end % of the search.  Return the best boundary and the associated MI & data