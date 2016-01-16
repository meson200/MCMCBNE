function [data_d,data_d_missing,mi,obj_disc] = kyu_discretize(data_c,data_c_missing,class,disc,bins,varargin)

% discretize data
% inputs
% data_c: continuous data 
% data_c_missing: continuous data with NaN
% disc: discretization method
% -disc = 1: median
% -disc = 2: supervised (maximum mutual information)
% -disc = 3: unsupervised (k-means)
% -disc = 4: apply the binning used for training set (parameters
% given in varargin)
% varargin: (for testing set only) a struct for discretization params created from training 

% outputs
% data_d: discretized data, missing values filled
% data_d_missing: discretized data, missing values NaN
% boundary_filled: bin boundaries for imputed data
% boundary_missing: bin boundaries for the raw data, NaN instances excluded
% mi: mutual information with a class (supervised discretization only)
% obj_disc: discretization field of the post processing object
% (obj_pp.disc)
% saves the following fields for later use:
% .disc = discretization method
% .bins = number of bins
% .boundary_filled: bin boundary(s), from lowest to highest, for imputed data
% .boundary_missing: bin boundary(s) for non-imputed data
% (boundary values are in normalization scale, so data needs to be
% normalized using obj_pp.norm before these boundaries are applied)


if ~isempty(varargin) 
    obj_disc = varargin{1,1}; 
else
    obj_disc = struct();
    obj_disc.disc = disc;
    obj_disc.bins = bins;
end


data_d = zeros(size(data_c));
data_d_missing = nan(size(data_c_missing));
boundary_filled = zeros(size(data_c,2),1);
boundary_missing = zeros(size(data_c,2),1);
mi = zeros(size(data_c,2),1);

for i=1:size(data_c,2)
    filled = data_c(:,i);
    unfilled = data_c_missing(:,i);
    nanindex = find(isnan(unfilled));
    pindex = find(~isnan(unfilled));
    
    if numel(unique(filled)) <= 2 % if the variable is already binary
        binned = filled;
        if min(binned) < 0
            binned(binned>0) = 2;
            binned(binned<0) = 1;
        end
        boundary_filled(i) = mean(filled);
        binned_missing = unfilled;
        if min(binned_missing) < 0
            binned_missing(binned_missing>0) = 2;
            binned_missing(binned_missing<0) = 1;
        end
        boundary_missing(i) = nanmean(unfilled);
        try
            mi(i) = MIarray(binned,class);
        catch
            mi(i) = NaN;
        end
    else
        if disc == 1 % median discretization
            boundary_filled(i) = median(filled);
            boundary_missing(i) = nanmedian(unfilled);
            binned = (filled >= boundary_filled(i))+1;
            binned_missing = (filled > boundary_missing(i))+1;
        elseif disc == 2 % MI discretization
            if bins == 3
                [~,~,binned,~] = opt3bin(filled,class); 
                [~,~,binned_missing,mi(i)] = opt3bin(unfilled(pindex),class(pindex));
            else
                binned_missing = nan(size(filled));
                [~,boundary_filled(i), binned] = opt2bin_kyu(filled,class,300);   
                [mi(i),boundary_missing(i), binned_temp] = opt2bin_kyu(unfilled(pindex),class(pindex),300);   
                binned_missing(pindex) = binned_temp;
            end 
        elseif disc == 3 % KM discretization
            sorted = sort(filled);
            %[binned,~] = kmeans(filled,2); %conventional
            [binned,~] = kmeans(filled,2,'Start',[sorted(1) sorted(length(sorted))]','EmptyAction','drop');
            %boundary_filled(i) = max(filled(binned==1));
            boundary_filled(i) = min(filled(binned==2));
            sorted_m = sort(unfilled);
            %[binned_missing,~] = kmeans(unfilled,2); % conventional
            [binned_missing,~] = kmeans(unfilled,2,'Start',[sorted_m(1) sorted(length(sorted_m))]','EmptyAction','drop');
            boundary_missing(i) = max(unfilled(binned_missing==1));
            try
                mi(i) = MIarray(binned,class);
            catch
                mi(i) = NaN;
            end
        elseif disc == 4 % apply bin boundaries in obj_disc
            binned = zeros(size(filled));
            for u = 1:numel(filled)
                binned(u) = (filled(u) > obj_disc.boundary_filled(i))+1;
                if ~isnan(unfilled(u)) 
                    binned_missing(u) = (unfilled(u) > obj_disc.boundary_missing(i))+1;
                else
                    binned_missing(u) = NaN;
                end
            end
        else    
            
        end
    end
    data_d(:,i) = binned;
    data_d_missing(:,i) = binned_missing;
    %[r,p(i)] = corr(unfilled(pindex),class(pindex));
end

obj_disc.boundary_filled = boundary_filled;
obj_disc.boundary_missing = boundary_missing;