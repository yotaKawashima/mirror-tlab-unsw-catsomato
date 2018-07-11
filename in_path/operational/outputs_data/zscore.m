function data_out = zscore(data)

% zscore: Finds the zscore of each value in a 3D matrix
%    data_out = zscore(data) finds the zscore using z = (x-mu)/std, with
%    mu and std calculated across each trial and channel.

% Written by Rannee Li, Dec 2014
% Last updated Jul 2015 (help file)

data_out = data;  % other fields are copied.

% preallocate
zscored = zeros(size(data.trial));

% calculate
nSamples = data.custom.nsamples;
for trial = 1:data.custom.ntrials
    
    % choose data
    x_sc = data.trial(:, :, trial); %ch
    
    % calculate
    mean_x = mean(x_sc, 2);
    mean_x_long = repmat(mean_x, 1, nSamples);
    
    std_x = std(x_sc,1,2); %, 1, nSamples);
    std_x_long = repmat(std_x, 1, nSamples);
    
    % find zscores
    sing_tr = (x_sc - mean_x_long) ./ std_x_long;
    
    % store
    zscored(:, :, trial) = sing_tr; %ch
    drawnow;
end


% rename
data_out.trial = zscored;

% change parameters
data_out.custom.filename = [data.custom.filename '_zscored'];
data_out.custom.datatype.iszscored = true;