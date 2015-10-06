function category = kyu_BN_RP_CategorizeVariables(labels)

% categorize variables into 6 categories as specified below:
% 1: physical 
% 2: biological-baseline + PTV volume
% 3: biological-post radiation
% 4: clincal factors 
% 5: RP

nvar = numel(labels);
category = [];
list_phy = {'MLD','MHD','V5','V10','V20','V30','Vhot','FracSize','MLD_BED'};
list_bio = {'a2m','TGF','IL','OPN','ACE','TNF'};
list_cli = {'age','smoking','IP','COPD','central','PTVCOMSI'};
list_spc = {'PTVvol'};

for i=1:nvar
    if occurrence(labels{i},list_phy) 
        cate = 1;
    elseif occurrence(labels{i},list_bio) 
        [root,suffix] = strtok(labels{i},'_');
        if strcmpi(suffix,'_pre')==1 
            cate = 2;
        else
            cate = 3;
        end
    elseif strcmpi(labels{i},'RP')~=0
        cate = 5;
    elseif strcmpi(labels{i},'PTVvol')~=0    
        cate = 2;
    else    
        cate = 4;
    end    
    category = [category cate];    
end