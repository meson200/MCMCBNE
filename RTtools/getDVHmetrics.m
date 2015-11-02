function [MLD,Vx,Vhot,gEUD] = getDVHmetrics(DVH,Dcutoff,Dhot,NumFrac,FracSize,DVHsettings,a)

% compute threeDVH metrics from DVH
%
% inptus:
% DVH: 'dvh' component of DVHobject (see kyu_readMIMDVH.m) 
%      column 1: dose (Gy) column 2: differential DVH (percentage vol.)
% Dcutoff: cutoff dose for Vx (low dose spillage)
% Dhot: cut off for Vhot (high dose spillage or hot spot)
% NumFrac/FracSize: number of fractions/fraction size
% DVHsettings: 
% a: weighting coefficient for the gEUD model

% outputs:
% MLD: mean dose to the structure
% Vx: percentage volume that received > Dcutoff
% Vhot: percentage volume that received > Vhot


if ~isempty(DVH)
    % LQ model based correction for fraction size
    if DVHsettings.NTDuse == 2
        DVH = DVH_normalize_NTD(DVH,NumFrac,FracSize,1,1,DVHsettings);
    end
    DVH_dose = DVH(:,1);
    DoseBinSize = DVH_dose(2)-DVH_dose(1);
    % dose values correspond to the left boundaries of bins
    % shift it to the central value of the bins
    DVH_dose_2 = DVH_dose(2:end);
    DVH_dose_2 = [DVH_dose_2; DVH_dose_2(end)];
    DVH_dose = (DVH_dose + DVH_dose_2)/2;
    % convert to daily dose
    %DVH_dose = DVH_dose/(NumFrac);
    
    DVH_volume = DVH(:,2);
    totalvolume = sum(DVH_volume);
    MLD = DVH_dose'*DVH_volume/totalvolume;
    gEUD = DVH_volume'*(DVH_dose.^a)/sum(DVH_volume);
    gEUD = gEUD^(1/a);
    [~,ind_Vx] = min(abs(DVH_dose-Dcutoff));
    [~,ind_Vhot] = min(abs(DVH_dose-Dhot));
    for i=1:length(DVH_volume)
       DVH_volume_cum(i) = (totalvolume - sum(DVH_volume(1:i)))/totalvolume;
    end
    % convert Vx and Vhot from absolute to relative volume
    Vx = DVH_volume_cum(ind_Vx)*100; 
    Vhot = DVH_volume_cum(ind_Vhot)*100;
    
else
    MLD = NaN;
    Vx = NaN;
    Vhot = NaN;
       
end