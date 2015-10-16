function DVH_out = DVH_normalize_NTD(DVH_in,NumFrac,FracSize,txnorm,clip,DVHsettings)

% Normalize x-axis (dose bins) of a DVH to 2-Gy-per fraction equivalent
% (EQD) using a LQ model 
% inputs:
% DVH_in: input DVH. column1: dose, column2: percentage dose
% FracSize/NumFrac: fraction size/num. of fractions
% txnorm: an option (on:2) to normalize dose to % of prescription dose 
% clip: to remove bins with zero volume (2)
% DVHsettings: 

abr = DVHsettings.abr;
Dhotmin = DVHsettings.HotspotDef;

DVH_dose = DVH_in(:,1);
DVH_volume = DVH_in(:,2);
Tx = FracSize*NumFrac*Dhotmin;
if DVHsettings.NTDuse == 2
    Tx_corr = Tx*(abr+Tx/NumFrac)/(abr+2);
else
    Tx_corr = Tx;
end
for i = 1:numel(DVH_dose)
    if DVHsettings.NTDuse == 2
        DVH_dose(i) = DVH_dose(i)*(abr+DVH_dose(i)/NumFrac)/(abr+2); % EQD2_4
    end
    if txnorm == 2
        DVH_dose(i) = DVH_dose(i)/Tx_corr*100;
    end
end
if clip == 2
    zeroindx = find(DVH_volume<10^-7);
    zeroindx = zeroindx(zeroindx>200);
    zeroindx = min(zeroindx);
    DVH_out = [DVH_dose(1:zeroindx-1) DVH_volume(1:zeroindx-1)];
else
    DVH_out = [DVH_dose DVH_volume];
end