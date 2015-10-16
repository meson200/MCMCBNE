function DVHobject = kyu_readMIMDVH(filepath,varargin)

% reads a dvh file written in a MIMvist agenerated .csv format
% and outputs the DIFFERENTIAL DVH
% (y axis: percentage volume)
% inputs:
% filepath: full path of the dvh file to read
% outputs: a list of DVH objects for all the structures containing the name
% in an input 'filter'
% varargin (optional) : a string filter applied to structure names.
% if varargin == 'LUNG', it imports only the structures which name contains
% 'LUNG'
%
% DVHobjects has the following fields:
% name: structure name 
% dvh: DVH array. 
%      column 1: dose (Gy) column 2: differential DVH (percentage vol.)
% volume: structure volume in cc

filter = '';
if ~isempty(varargin) 
    filter = varargin{1,1}; 
end

fid = fopen(filepath,'r');
% scan the header file
% skip the first line
line = fgetl(fid);
line = fgetl(fid);
header = textscan(line, '%s','delimiter',',');
header = header{1}(2:end);
numstr = numel(header);
textformat = repmat('%s',1,numstr+1); 

% read the numeric part
num = textscan(fid, textformat, 'delimiter',',');

cumulative = false; % convert to differential!!!
% find the relevant columns
% find the columns which name includes a string in the 2nd input arg.
col = [];
strnames = [];
strvol = [];
for i=1:numstr
    tmpstr = header{i};
    tmpstr2 = textscan(tmpstr,'%s%s%s','delimiter',' ');
    tmpname= tmpstr2{1}{1};
    tmpvol = textscan(tmpstr2{3}{1},'%s%s','delimiter',')');
    tmpvol = tmpvol{1}{1};
    tmpname(regexpi(tmpname,'"'))  = [];
    str_name = strtrim(tmpname);
    aaa = regexpi(str_name,filter);
    if isempty(filter) || ~isempty(aaa)
       col = [col i]; 
       strnames = [strnames cellstr(str_name)]; 
       strvol = [strvol str2num(tmpvol)];
    end
end

DVHobject = struct('name',{'a'},'dvh',[1],'volume',[1]);

for p = 1:length(col)
  
    DVH_dose = num{1};
    DVH_volume = num{col(p)+1};
    try
        DVH_dose = convert_to_num(DVH_dose);
    catch exception 
        getReport(exception)
    end
    DVH_volume = convert_to_num(DVH_volume);

    if cumulative 
        DVH_volume_2 = DVH_volume;
    else % convert cumulative to differential 
        DVH_volume_2 = zeros(size(DVH_volume));
        DVH_volume = [DVH_volume; 0];
        for i = 1:length(DVH_dose)
            DVH_volume_2(i) = DVH_volume(i)-DVH_volume(i+1);
        end
    end
    dvhtemp = [DVH_dose DVH_volume_2];
    dvh{p} = dvhtemp;
        
end
 
DVHobject.name = strnames;
DVHobject.dvh = dvh;
DVHobject.volume = strvol;

%DVHobject = struct('name',strnames,'dvh',dvh);

fclose(fid);

function num = convert_to_num(str)

counts = length(str);
num = zeros(counts,1);
for i = 1:counts
    this = str(i);
    this = this{1};
    if isempty(this)
        num(i) = 0;
    else
        if this(1) == '"'
            if length(this)>2
                this = this(2:end-1);
                num(i) = str2num(this);
            else
                num(i) = 0;
            end
        else    
            num(i) = str2num(this);
        end
    end
end
        
    
    
    

