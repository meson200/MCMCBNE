function output = kyu_estimateNTCP(studyID,FracSize,COMSI,source)

% calculate dosimetric variables & RP risks from DVH
% different formats and directories are used for DVHs from two institutions
% inputs:
% studyID: from kyu_BN_readdata_combined 
% FracSize: from kyu_BN_readdata_combined
% source: (1) WashU (2) McGill 
% output: [MLD V20 gEUD Bradley]

% LQ model based correction for fraction size
abr = 4; % alpha-beta ratio for pneumonitis (Bentzen et al. 2007)
fac = (FracSize+abr)/(2+abr);


if studyID < 10
    filename = cat(2,'L00',num2str(studyID),'_dvh');
else
    filename = cat(2,'L0',num2str(studyID),'_dvh');
end

if source==1
    if studyID == 7
       filename = cat(2,filename,'_cum'); 
    end
    filename = cat(2,filename,'.xls');
    directoryname = '/Users/kyu/Desktop/Patients/NSCLC_WUSTL/DVH/';
    fullpath = cat(2,directoryname,filename);
   
else
    directoryname = '/Users/kyu/Desktop/Patients/NSCLC_McGill/DVH/';
    fullpath = cat(2,directoryname,filename);
end

filethere = exist(fullpath,'file');
if filethere
    
    if source == 1
        if studyID == 7
            [DVH_dose,DVH_volume] = kyu_readEclipseDVH_WashU(filename,2);
        else
            [DVH_dose,DVH_volume] = kyu_readEclipseDVH_WashU(filename,1);
        end
    else
        [DVH_dose,DVH_volume] = kyu_readEclipseDVH_McGill(filename);
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
