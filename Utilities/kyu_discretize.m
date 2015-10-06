function [data_d,data_d_missing,boundary_filled,boundary_missing,mi] = kyu_discretize(data_c,data_c_missing,class,disc,bins)

% discretize data
% data_d: discretized data, missing values filled
% data_d_missing: discretized data, missing values NaN
% boundary_filled: bin boundaries for imputed data
% boundary_missing: bin boundaries for the raw data, NaN instances excluded
% mi: mutual information with a class (supervised discretization only)


data_d = zeros(size(data_c));
data_d_missing = zeros(size(data_c_missing));
boundary_filled = zeros(size(data_c,2),1);
boundary_missing = zeros(size(data_c,2),1);

for i=1:size(data_c,2)
    filled = data_c(:,i);
    unfilled = data_c_missing(:,i);
    nanindex = find(isnan(unfilled));
    pindex = find(~isnan(unfilled));
    if disc == 2
        if bins == 3
            [~,~,binned,~] = opt3bin(filled,class); 
            [~,~,binned_missing,mi(i)] = opt3bin(unfilled(pindex),class(pindex));
         else
            [~,boundary_filled(i), binned] = opt2bin(filled,class,300);   
            [mi(i),boundary_missing(i), binned_missing] = opt2bin(unfilled(pindex),class(pindex),300);   
        end 
        data_d_missing(pindex,i) = binned_missing;
        data_d_missing(nanindex,i) = NaN*ones(1,numel(nanindex));
    else
        if numel(unique(filled)) == 2 % if the variable is already binary
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
        else
            sorted = sort(filled);
            %[binned,~] = kmeans(filled,2); %conventional
            [binned,~] = kmeans(filled,2,'Start',[sorted(1) sorted(length(sorted))]','EmptyAction','drop');
            boundary_filled(i) = max(filled(binned==1));
            sorted_m = sort(unfilled);
            %[binned_missing,~] = kmeans(unfilled,2); % conventional
            [binned_missing,~] = kmeans(unfilled,2,'Start',[sorted(1) sorted(length(sorted))]','EmptyAction','drop');
            boundary_missing(i) = max(unfilled(binned_missing==1));
        end
        mi(i) = MIarray(binned,class);
        data_d_missing(:,i) = binned_missing;
    end
    data_d(:,i) = binned;
    [r,p(i)] = corr(unfilled(pindex),class(pindex));
end