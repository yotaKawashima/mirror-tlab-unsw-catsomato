function data_out = dwnsmpl(data, target)

% dwnsmpl: downsamples the data
% 
%   data_out = dwnsmpl(data, target) takes data in struct and resamples it
%   to the timepoints that are in target. It chooses the closest points
%   in data.time to target. data_out is a truct with the same fields as
%   data.

% Written my Rannee Li, Jan 2015
% Last updated Jul 2015 (data.time updated to 1x1 cell)

data_out = data;

rep_target = repmat(target, data.custom.nsamples, 1);
rep_time = repmat(data.time{1}', 1, length(target));

[~, inds] = min(abs(rep_time - rep_target));
% inds = unique(inds);

% for i = 1:data.custom.ntrials
%     data_out.time{i} = target;
% end

data_out.time{1} = target;

data_out.trial = data.trial(:, inds, :);


% fixing other variables
data_out.fsample                    = 1/(target(2) - target(1));

data_out.custom.nsamples            = length(target);
data_out.custom.filename            = [data.custom.filename '_dwnsmpl'];
data_out.custom.datatype.isdwnsmpl  = true;
data_out.custom.true_time           = data.time{1}(inds);