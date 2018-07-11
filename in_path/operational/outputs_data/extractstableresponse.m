function data_out = extractstableresponse(data, timepoints)
% extractstableresponse: selects a portion of data
%
%   data_out = extractstableresponse(data, timepoints) extracts from
%   data.trial only the timepoints that lie between the first and last
%   values of timepoints
%   
%   The function detects if frequency or time data is entered.

% Written by Rannee Li, Dec 2014
% Last updated Jul 2015 (fixed time cell size)
% To update: error and datatype checking may be too complex

data_out = data;  % other fields are copied.

% error checking
if numel(timepoints)<=1 || numel(timepoints)~=floor(numel(timepoints))
    error('numel(timepoints) must be a positive integer greater than 1.')
end

[~, c, ~] = size(data.trial);

iffreq = false;
try
    timestamp = data.time{1};
catch 
    warning('No times detected. Trying frequencies instead.')
    timestamp = data.freq{1};
    iffreq = true;
end

if length(timestamp)~=c
    error('length(timestamp)~=size(data, 2)')
end

iscellarr = istype(data.custom.datatype, 'cellarr');
if iscellarr
    error('Data in cell array structure and cannot be processed.')
end


% find the element in timescale that is the closest to timepoints(1)
find_tmp = abs(timestamp - timepoints(1));
ind(1) = find(find_tmp==min(find_tmp));

% find the same for timepoints(end)
find_tmp = abs(timestamp - timepoints(end));
ind(2) = find(find_tmp==min(find_tmp));

data_out.custom.nsamples = ind(2)-ind(1)+1;

timestamp_new = timestamp(ind(1):ind(2));

% extract
data_out.trial = data.trial(:, ind(1):ind(2), :);




if iffreq
    % reallocate data.time
%     data_out.freq = cell(1, data.custom.ntrials);
%     for i = 1:data.custom.ntrials
%         data_out.freq{i} = timestamp_new;
%     end
    data_out.freq{1} = timestamp_new;

    % set type
    data_out.custom.datatype.isfxtoxxx = true;
    data_out.custom.filename = [data.custom.filename '_f' num2str(timepoints(1)) 'to' num2str(timepoints(2))];
    
else    
    % reallocate data.time
%     data_out.time = cell(1, data.custom.ntrials);
%     for i = 1:data.custom.ntrials
%         data_out.time{i} = timestamp_new;
%     end
    data_out.time{1} = timestamp_new;
    
    % set type
    data_out.custom.datatype.isevkresp = true;
    data_out.custom.filename = [data.custom.filename '_evkresp'];
end