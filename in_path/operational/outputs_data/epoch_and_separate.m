function data = epoch_and_separate(filename_in, tres, timestamp, lfp1_data, lfp2_data, stimFreqY, stimFreqX, stimMeshY, stimMeshX, stimTime, output_dest, filename_out) %#ok - used in eval line 43

% epoch_and_separate: Epochs the data from into 3D matrix.
%    epoch_and_separate(...) epochs the data (for this format) into a
%    struct, data, and saves it to the folder epoched. For more information
%    on the output see the supporting documentation.
%
%    data = epoch_and_separate(...) returns the output structure to the
%    workspace.
%   
%    The argument list: filename_in, tres, timestamp, lfp1_data, lfp2_data, 
%       stimFreqY, stimFreqX, stimMeshY, stimMeshX, stimTime, output_dest
%       
%    Due to the different types of filenames, an additional argument,
%    filename_out can be added to the end of the list to define the start 
%    of the output filename.
%
%    Note that the file must be loaded BEFORE this function is run.

% REMEMBER TO PRELOAD DATA

% Written by Rannee Li, Dec 2014. 
% Last updated Jul 2015.

fprintf('\nEntered function epoch_and_separate.\n')

%% write variables to the data structure
% 

datatype.iszscored = false;
datatype.isbiprref = false;
datatype.ispwrspec = false;

data.fsample                = 1/(str2double(num2str(tres)));        % 1x1 double
data.custom.datatype        = datatype;                             % struct

%% deal with the filename
if nargin < 11
catnum = filename_in(4:11);
recnum = filename_in(23);
ttype  = filename_in(28:31);
CatName = ['C' catnum '_R' num2str(str2double(recnum), '%02i') '_T' ttype];
else 
    CatName = filename_out;
end
%% epoch data

fprintf('Epoching data.\n')


for S = [1, 2];
    
    fprintf('\tEpoching S%i\n', S)
    
    % rename the lfp
    eval(['lfp = lfp' num2str(S) '_data;']);
    
    % find nConds
    [CondHeight, CondWidth] = size(lfp); %#ok<USENS> % warning in this line resolved by the eval in line 39
    data.custom.subplotconfig   = [CondHeight, CondWidth];          % 1x2 double
    nConds = CondHeight * CondWidth;
    data.custom.conditions      = [ 0 nConds];                      % 1x2 double
    
    % find nChannels
    
    [ConfigRows, ConfigColumns, nSamples] = size(lfp{1}{1});
    nChannels = ConfigRows * ConfigColumns;
    
    data.custom.nsignals        = nChannels;                        % 1x1 double
    data.custom.nsamples        = nSamples;                         % 1x1 double
    data.custom.spatialconfig   = [ConfigRows, ConfigColumns];      % 1x2 double
    data.custom.area            = S;                                % 1x1 double
    
    % pre-allocate
    CondNames = cell(16, 1);
    chlabels = cell(1, data.custom.nsignals);
    
    AreaName = [CatName '_S' num2str(S) '_'];
    
    for cond = 1:nConds
        % find name
        CondNames{cond} = ['F' num2str(stimFreqY, '%03i') 'A' num2str(stimMeshY(cond), '%03i') ...
            '_F' num2str(stimFreqX, '%03i') 'A' num2str(stimMeshX(cond), '%03i')];
        
        % take out the unepoched matrix
        unepoch = lfp{cond};
        
        [nTrials, ~] = size(lfp{cond});
        data.custom.ntrials         = nTrials;      % 1x1 double
        
        % allocate
        epoched = zeros(nChannels, nSamples, nTrials);
       
        % epoch
        ch = 1;
        for electrode_x = 1:ConfigRows
            for electrode_y = 1:ConfigColumns
                
                for trial = 1:nTrials
                    if ch<=nChannels
                        epoched(ch, :, trial) = unepoch{trial}(electrode_y, electrode_x, :);
                        
                        % create chlabels
                        if trial == 1
                            chlabels{ch} = ['rawdat_ch' num2str(ch, '%03i')]; 
                        end
                    else
                        error('channel number invalid')
                    end
                end
                ch = ch+1;
            end
        end
        fileName = [AreaName CondNames{cond} '_epoched'];

        data.custom.filename        = fileName;                     % string
        
        % create time
        stimStart = stimTime.rampup + stimTime.presine; % so that t=0 is the onset of vibration
        data.time{1} = timestamp-stimStart;  % EDIT 06/07: changed this to a 1x1 cell. See documentation
        
        data.label = chlabels;                                      % 1xnChannels cell (string)
        data.custom.conditions(1) = cond;
        
        % save
        data.trial = epoched;
        eval(['save ' fileName ' data'])
        
        
    end
    
    % move everything
    eas_filemover(AreaName, output_dest)
    

end

if nargout < 1
    clear data
end

fprintf('Epoching data complete.\n')
fprintf('Exited function epoch_and_separate.\n')

function eas_filemover(files_to_move, output_dest)
% adapted from afpc_filemover

% deal with spaces in file path
k = strfind(output_dest, ' ');

if ~isempty(k)
    warning('Spaces detected in path name')
    for a = 1:numel(k)
        output_dest(k(a):end+1) = ['\' output_dest(k(a):end)];
    end
end

% find out if the folder exists.
if exist(output_dest, 'dir')~=7
    % directory does not exist, so make it. 
    warning('Destination directory not found. Making directory.')
    system(['mkdir ' output_dest]);
end

% move the files
system(['mv ' files_to_move '* ' output_dest]);