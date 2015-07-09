function output = kyu_estimateNTCP(studyID,FracSize,COMSI)

% calculate dosimetric variables & RP risks from DVH
% different formats and directories are used for DVHs from two institutions
% inputs:
% studyID: from kyu_BN_readdata_combined 
% FracSize: from kyu_BN_readdata_combined
% source: (1) WashU (2) McGill (3) CHUM
% output: [MLD V20 gEUD Bradley]

% LQ model based correction for fraction size
abr = 4; % alpha-beta ratio for pneumonitis (Bentzen et al. 2007)
fac = (FracSize+abr)/(2+abr);

dirname = '~/Box Sync/SKyu/lungdata/DVH/';
cd(dirname);
filename = strcat(studyID,'_dvh');
inst = studyID(1);
switch inst
    case 'W' % WashU
        ext = '.xls';
    case 'L'
        ext = '';
    case 'C'
        ext = '.csv';
end
fullpath = strcat(2,dirname,filename,ext);
filethere = exist(fullpath,'file');
if filethere
    
    switch inst
    case 'W' 
        if studyID == 7
            [DVH_dose,DVH_volume] = kyu_readEclipseDVH_WashU(fullpath,2);
        else
            [DVH_dose,DVH_volume] = kyu_readEclipseDVH_WashU(fullpath,1);
        end
    case 'L'
        [DVH_dose,DVH_volume] = kyu_readEclipseDVH_McGill(fullpath);
    case 'C'
        [DVH_dose,DVH_volume] = kyu_readEclipseDVH_CHUM(fullpath);
    end

    % calculate MLD first
    totalvolume = sum(DVH_volume);
    MLD = DVH_dose'*DVH_volume/totalvolume;
    %disp(MLD)

    % calculate V20 and V30
    for i=1:length(DVH_volume)
       DVH_volume_cum(i) = (totalvolume - sum(DVH_volume(1:i)))/totalvolume;
    end
    
    [~,ind_V20] = min(abs(DVH_dose-20/fac));
    [~,ind_V30] = min(abs(DVH_dose-30/fac));
    V20 = DVH_volume_cum(ind_V20)*100;
    V30 = DVH_volume_cum(ind_V30)*100;

    % gEUD
    % used n=1.03 (quantec)

    n = 1.03;
    gEUD = 0;
    gEUD_NTD = 0;
    for i=1:numel(DVH_dose)
        gEUD = gEUD + DVH_volume(i)*(DVH_dose(i)*fac)^(1/n);
    end 
    gEUD = (gEUD/totalvolume)^n;

    % NTCP_LKB
    % used Burman's parameter (1991)
    TD50 = 24.5;
    m = 0.18;
    t = (gEUD - TD50)/(m*TD50);
    NTCP_LKB = normcdf(t,0,1);

    % Bradley et al. (2007)
    eta = [1 MLD COMSI] * [-1.5 0.11 -2.8]';
    Bradley = drxlr_invlogit(eta);
        
    disp([num2str(studyID),' MLD:',num2str(MLD*100),' V20:',num2str(V20),' V30:',num2str(V30)]);

else
    MLD = NaN;
    V20 = NaN;
    V30 = NaN;
    gEUD = NaN;
    Bradley = NaN;
end
output = [MLD V20 gEUD Bradley];
