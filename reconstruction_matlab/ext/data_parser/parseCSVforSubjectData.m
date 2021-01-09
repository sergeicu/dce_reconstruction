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
function [  files, timePoints, patientName , MRN, refNii, badCoils ] = parseCSVforSubjectData( patientName, csvFile )

    %% read subjects list file containing patient name/date and address and name of .dat files 
    
    %% read raw CSV file
    raw_csv = textread(csvFile,'%s','delimiter',',','emptyvalue',NaN);

    % find end of arrays with the "END" marker
    eof = strfind(raw_csv,'END');
    end_ix = find(not( cellfun('isempty', eof) ));


    % reshape cell array. Use only the columns/rows that are consistent
    coherent_rows = min(find(mod(end_ix,end_ix(1)) >0 ))-1;      % if each row modulo first is 0 we are good. Otherwise, there's a mess
    if numel(coherent_rows) == 0
       coherent_rows = numel(raw_csv) / end_ix(1); 
    end
    table_csv = reshape( raw_csv(1:coherent_rows*end_ix(1)) ,  [end_ix(1), coherent_rows] );

    % select subject names
    [startNames_ix]=find(~cellfun(@isempty,strfind(table_csv(:,1),'Names'))) +1;
    names_list = table_csv(startNames_ix:end-1,1);
    fprintf('Names position is in the entry %d. This should be a 2.\n',startNames_ix);

    % if no patient specified, take them all!
    if numel(patientName) == 0
        patientName = names_list;
    end


    %% get pointers to columns
    % retrieve columns with file info
    [~,file_info_cols] = find(not( cellfun('isempty', strfind(table_csv,'reconFile') ) ));
    % retrieve columns with time points info
    [~,time_pts_cols] = find(not( cellfun('isempty', strfind(table_csv,'timePoints') ) ));
    % retrieve columns with info about MRNs
    [~,MRN_cols] = find(not( cellfun('isempty', strfind(table_csv,'MRN') ) ));
    
    
    % retrieve columns with info about MRNs
    [~,refNii_cols] = find(not( cellfun('isempty', strfind(table_csv,'referenceNii') ) ));
    % retrieve columns with info about MRNs
    [~,badCoils_cols] = find(not( cellfun('isempty', strfind(table_csv,'bad coils') ) ));

    files = cell(1, numel(patientName));
    timePoints = cell(1,numel(patientName));
    MRN = cell(1,numel(patientName));
    refNii = cell(1,numel(patientName));
    badCoils = cell(1,numel(patientName));
    % for every patient requested
    for pt = 1:numel(patientName)

        % retrieve patient from list
        [row_number] = find(~cellfun(@isempty,strfind(names_list,patientName{pt}))) +2;      % adding 2 to be consistent with table (there is a '' and '"Names"' at the beginning)

        % GET .DAT FILES
        files_ix = ~cellfun( @isempty, table_csv( row_number, file_info_cols ) );
        files{pt} = table_csv( row_number, file_info_cols(files_ix) )';
        numOfDatFiles = numel(files{pt});%length(find(~cellfun( @isempty, table_csv( row_number, file_info_cols ) )));        
        
%         % clean up the " thingys
%         for ff =1:numOfDatFiles
%            files{pt}{ff} = char(extractBefore(extractAfter(files{pt}{ff},'"'),'"'));
% %            files{pt}{ff} = sprintf('%s',files{pt}{ff});
%         end
        
        
        % GET TIME POINTS
        time_pts_ix = ~cellfun( @isempty, table_csv( row_number, time_pts_cols ) );
        timePoints{pt} = cellfun( @str2num, (  table_csv( row_number, time_pts_cols(time_pts_ix)  ) ) );

        % GET MRN
        MRN{pt} = str2num(cell2mat( table_csv( row_number, MRN_cols )  ));

        % get referenceNii
        if numel(refNii_cols) > 0
            tmp =  table_csv( row_number, refNii_cols );
            refNii{pt} = tmp{1};
        end

        % get bad coils
        if numel(badCoils_cols) > 0
            temp = table_csv( row_number, badCoils_cols );
            eval(sprintf('badCoils{%d} = [%s];',pt, temp{1}));
        end
    end

end