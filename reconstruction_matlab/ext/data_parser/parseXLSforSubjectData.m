%%
%       This function parses the XLS file with all th subjects and returns the path to the desired files and the time points.
%
%       INPUTS:
%           - patientName - cell - containing the name of each subject
%                                       - str - name of the subject
%           -   xlsFile - str - path to the XLS file containing all the patients information
%
%       OUTPUTS:
%           - files - cell - each cell represents a subject and contains a cell with the paths to the data files
%                                   - cell - containing str to data path
%           - timePoints - cell - each cell represents a subject
%                                   - <1,NF>int - each entry represents the number of volumes in each data file.
%
%
%
function [  files, timePoints, patientName , MRN ] = parseXLSforSubjectData( patientName, xlsFile )

    %% read subjects list file containing patient name/date and address and name of .dat files 
    [~,txt,raw]  = xlsread(xlsFile);

    % if no patient specified, take them all!
    if numel(patientName) == 0
        patientName = txt(1,3:end);
    end

    files = cell(1, numel(patientName));
    timePoints = cell(1,numel(patientName));
    MRN = cell(1,numel(patientName));
    for pt = 1:numel(patientName)

        [rN,cN]=find(~cellfun(@isempty,strfind(txt,patientName{pt})));
        cN
        txt(:,cN)
        raw(:,cN)
        numOfDatFiles=length(find(~cellfun(@isempty,txt(4:8,cN))));

        files{pt}=txt(4:4+numOfDatFiles-1,cN);                       % retrieve .dat files with K-space data (possibly more than one)
        timePoints{pt} = cell2mat(raw(8:8++numOfDatFiles-1,cN));    % retrieve the number of time points / volumes to use

        MRN{pt} = cell2mat(raw(3,cN));
    end

end