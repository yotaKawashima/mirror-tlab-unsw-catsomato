function data_out = chronux_pwrspec_v2(data, params)

% chronux_pwrspec: Finds the power spectrum across each trial and channel
%   data_out = chronux_pwrspec(data) uses the mtspectrumc function (in the
%   chronux suite) to compute the power spectrum using default settings. It
%   calculates each trial and channel individually.
%
%   data_out = chronux_pwrspec(data, params) uses params to configure
%   mtspectrumc instead
%
%   Version 2 uses mtspectrumc's ability to compute a number of trials
%   together for one channel. (Compute power for each trial.)

% Modified from pwrspec by Rannee Li, Jul 2015
% Last edited Jul 2016 by Rannee Li


fprintf('\nEntered function chronux_pwrspec.\n')

% copy everything
data_out = data;

% rename the variables
nChannels = data.custom.nsignals;
%nTrials = data.custom.ntrials;
fsample = data.fsample;

% set chronux parameters if none provided
if nargin < 2 || isempty(params)
    params.Fs = fsample;
    params.tapers = [1 1];
    params.pad = 1;
end

fprintf('Calculating power spectrum for %s.\n', data.custom.filename)
% calculate by iteration
for ch = 1:nChannels
    
        % chronux requires format to be samples x trials
        timetrial = permute(data.trial(ch, :, :), [2, 3, 1]);
        
        % calculate
        [p_part,f_part] = mtspectrumc(timetrial, params);
        
        % store
        power_tmp = 10*log10(p_part'); 
        powers(ch, :, :) = permute(power_tmp, [3, 2, 1]); %#ok<AGROW> don't know how many freqs
        
        if ch==1
            freq = f_part';
        end
        
    
end
fprintf('Calculating complete.\n')

% place into struct, fix some vars
data_out.trial = powers;
data_out.custom.filename = [data.custom.filename '_cmtspwr'];
data_out.custom.datatype.ispwrspec = true;
data_out.custom.time = data.time{1}; % this is needed for some of the anova analysis

% copy freqs
data_out = rmfield(data_out, 'time');


data_out.freq{1} = freq;


data_out.custom.params = params;


if nargout < 1
    clear data;
    data = data_out; %#ok<NASGU> used in save below
    save(data_out.custom.filename, 'data')
    clear data_out
end
fprintf('Exited function chronux_pwrspec.\n')
