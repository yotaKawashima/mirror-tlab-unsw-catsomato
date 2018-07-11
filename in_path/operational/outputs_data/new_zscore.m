function data_out = new_zscore(data)
% new_zscore: Finds the zscore of each value in a 3D matrix
%    data_out = new_zscore(data) finds the zscore using z = (x-mu)/std, 
%    with mu and std calculated across each channel.
% 
%    Note that this is different from zscore as it operates for all trials
%    in a single channel instead of across each single trial in a single
%    channel.

% Written by Rannee Li, Dec 2014
% Last updated Jul 2015 (help file)

data_out = data; % copy all variables

% find the mean across trials 
mean_trials = mean(data.trial, 3);

% find the std across trials
std_trials = std(data.trial, 1, 3);

% 
data_out.trial = (data.trial - repmat(mean_trials, [1, 1, data.custom.ntrials])) ./ repmat(std_trials, [1, 1, data.custom.ntrials]);

% change parameters
data_out.custom.filename = [data.custom.filename '_newzscr'];
data_out.custom.datatype.isnewzscr = true;


% tmp1 = [];
% tmp2 = [];
% 
% % first put all trials into a single matrix
% data_trial = data.trial; % check this with Lisandro
% for tr = 1 : length(data_trial) 
%         tmp1 = [tmp1 ; data_trial{tr}(1,:)]; %#ok
%         tmp2 = [tmp2 ; data_trial{tr}(2,:)]; %#ok
% end
%
% if normalize_channels
%     
%     % Normalized channels
%     if along_trials
%         channel1 = (tmp1 - repmat(mean(tmp1), [nTrials 1]) )  ./ repmat(std(tmp1), [nTrials 1]);
%         channel2 = (tmp2 - repmat(mean(tmp2), [nTrials 1]) )  ./ repmat(std(tmp2), [nTrials 1]);
%     else
%         channel1 = (tmp1 - repmat(mean(tmp1, 2), [1 size(tmp1, 2)]) )  ./ repmat(std(tmp1), [nTrials 1]); 
%         channel2 = (tmp2 - repmat(mean(tmp2, 2), [1 size(tmp2, 2)]) )  ./ repmat(std(tmp2), [nTrials 1]);
%     end
%     data.trial = [];
%     for tr = 1 : nTrials
%         data.trial{tr} = cat(1, channel1(tr, :), channel2(tr, :));
%     end
%     clear tmp1 tmp2
%     disp('Channels Normalized')
% end