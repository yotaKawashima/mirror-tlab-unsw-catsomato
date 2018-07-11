function data = anp_data_around100(data_dir, filename_header, frange, ifunipol,TargFreq)

% anp_data: arranges data for the anova2 pipeline
%   data = anp_data(data_dir, filename_header) produces data.y, which can
%   be fed into anp_analyse for analysis, based on the file specified by
%   filename_header in the directory data_dir. It selects the data for all
%   channels present and all frequencies.
%
%   data = anp_data(data_dir, filename_header, frange) runs only the
%   frequencies given by frange. If frange is a vector, then the first and
%   last elements are considered to be the bounding frequencies. If frange
%   is a number it is treated as the maximum frequency to be considered. By
%   default is is all of the calculated frequencies
%
%   data = anp_data(data_dir, filename_header, frange, ifunipol) runs only
%   unipolar channels if ifunipol is true, and only bipolar channels if it
%   is false. By default it is false.
%
%   This function loads and saves data within itself.

% Written by Rannee Li, July 2015, editted by Phoebe You 22 Feb 2016 for analysing
% interesting response from 50Hz to 150Hz 

%% set defaults
if nargin < 4 || isempty(ifunipol)
    ifunipol = false;
end

%% find the file location and the number of conditions to be loaded
fnames = dir(fullfile(data_dir, [filename_header, '*.mat']));
nConds = length(fnames);

%% load and process the first dataset
% this allows automatic reading of the number of trials and the format of
% the output data
load(fullfile(data_dir, fnames(1).name))

% find number and arrangement of conditions
nTrials = data.custom.ntrials; %#ok<NODEF> data is a loaded variable

% find out which channels to select (only bipolar referencing channels)
lastunipol = prod(data.custom.spatialconfig);
if ifunipol
    chan = 1:lastunipol;
    pollet = '1';
else
    chan = lastunipol+1:data.custom.nsignals;
    pollet = '2';
end
label = data.label(chan);
nchan = numel(label);

% % find the frequencies
% if numel(frange)==1
%     [~, fmax_ind] = find_closest(data.freq{1}, frange);
%     fmin_ind = 1;
%     frange(2) = frange;
%     frange(1) = 0;
% else
%     [~, fmax_ind] = find_closest(data.freq{1}, frange(end));
%     [~, fmin_ind] = find_closest(data.freq{1}, frange(1));
% end

%% frange = [10 5] --- [segment values, number of segment pairs] eg,get the mean of 90-110,80-120,70-130,60-140,50-150
% the segment value is 10, and the number of segment pairs is 5
for iSeg = 1:frange(end)
    RangeFreq(iSeg,:) = [TargFreq-frange(1)*iSeg,TargFreq+frange(1)*iSeg];
    [~, fmax_ind(iSeg)] = find_closest(data.freq{1}, RangeFreq(iSeg,end));
    [~, fmin_ind(iSeg)] = find_closest(data.freq{1}, RangeFreq(iSeg,1));

end

% preallocate y, the output
subpcfg = data.custom.subplotconfig;

% Get the mean of power spectrum for one trial, one channel across the
% sepecific frequency range

for c = 1:nchan
    for t = 1:nTrials
        for s = 1:size(RangeFreq,1)
        NewData.trial(c,s,t) = mean(data.trial(c,fmin_ind(s):fmax_ind(s),t));
        end
    end
end

dim1 = nTrials*subpcfg(1); % increasing 200 Hz, across trials
dim2 = subpcfg(2); % increasing 23 Hz
% dim3 = fmax_ind-fmin_ind+1; % frequencies 
dim3 = frange(end); % number of segments pairs
dim4 = nchan; % channels

y = zeros(dim1, dim2, dim3, dim4);

namematrix = cell(subpcfg); % the order the data in organised in 

% put variables in places!
k = 1;
col = mod(k-1, subpcfg(2)) +1; % column and row where the data goes
row = (floor((k-1)/subpcfg(1))*nTrials) +1; 

% dat = data.trial(chan, fmin_ind:fmax_ind, :);
dat = NewData.trial;

y(row:row+nTrials-1, col, :, :) = permute(dat, [3, 4, 2, 1]);

namematrix{floor((k-1)/subpcfg(1))+1, col} = data.custom.filename(24:40);


%% load and process rest of data
for k = 2:nConds
    %% load file
    load(fullfile(data_dir, fnames(k).name))
    clear NewData.trial nTrials
    %% extract only the bipolar data at given frequencies, save into the matrix
    nTrials = data.custom.ntrials; % some condition has 7 trials but some has 6 trials
    
    col = mod(k-1, subpcfg(2)) +1;
    row = (floor((k-1)/subpcfg(1))*nTrials) +1;
    
    TargFreq = 100;
    for iSeg = 1:frange(end)
        RangeFreq(iSeg,:) = [TargFreq-frange(1)*iSeg,TargFreq+frange(1)*iSeg];
        [~, fmax_ind(iSeg)] = find_closest(data.freq{1}, RangeFreq(iSeg,end));
        [~, fmin_ind(iSeg)] = find_closest(data.freq{1}, RangeFreq(iSeg,1));
        
    end
    % Get the mean of power spectrum for one trial, one channel across the
    % sepecific frequency range
    for c = 1:nchan
        for t = 1:nTrials
            for s = 1:size(RangeFreq,1)
                NewData.trial(c,s,t) = mean(data.trial(c,fmin_ind(s):fmax_ind(s),t));
            end
        end
    end
    
    %     dat = data.trial(chan, fmin_ind:fmax_ind, :);
    dat = NewData.trial;
    y(row:row+nTrials-1, col, :, :) = permute(dat, [3, 4, 2, 1]);
    
    namematrix{floor((k-1)/subpcfg(1))+1, col} = data.custom.filename(24:40);
end

%% collect variables, save

% variables in the 'main menu'
dat_tmp.conditions  = namematrix;
dat_tmp.chanlabel   = label;
dat_tmp.freq{1}     = data.freq{1}(fmin_ind(end):fmax_ind(end));

% data.custom updates
dat_tmp.custom                      = data.custom;
dat_tmp.custom.datatype.isanovadata = true;
dat_tmp.custom                      = rmfield(dat_tmp.custom, {'conditions', 'nsamples'});
dat_tmp.custom.unipolar             = ifunipol;
dat_tmp.custom.nsignals             = nchan;
dat_tmp.custom.filename(24:40)      = ['F' num2str(RangeFreq(frange(end),1), '%03i') 'to' ...
    num2str(RangeFreq(frange(end),2), '%03i') '_P' pollet '_anov'];
dat_tmp.custom.filename             = [dat_tmp.custom.filename '_adatain'];

% delete comm
data = dat_tmp;
data.y          = y;

% save
save(data.custom.filename, 'data')

if nargout < 1
    clear data
end
