function timecourse(data_dir, cat_name, area, winsize, params)


% Partially adapted from chronux_pwrspec_v2, 13 Nov 2017

fprintf('\nEntered function timecourse.\n')

% set defaults
if nargin < 3
    area = 'S2';
end
if nargin < 4
    winsize = 0.5;
end

dir_sec = fullfile(data_dir, cat_name, 'epoched_rsampsl_biprref/');

loadname = dir(fullfile(dir_sec, ['*' area '*.mat']));

for k = 1:numel(loadname);
    load(fullfile(dir_sec, loadname(k).name))
   
    clear powers
    
    % only care about the bipolar channels!!!
    if numel(data.label)>prod(data.custom.spatialconfig)
        c1 = prod(data.custom.spatialconfig)+1;
        data.label = data.label(c1:end);
        data.trial = data.trial(c1:end, :, :);
        data.custom.nsignals = numel(data.label);
    end
    
    % copy everything
    data_out = data;
    
    % rename the variables
    nChannels = data.custom.nsignals;
    %nTrials = data.custom.ntrials;
    fsample = data.fsample;
    
    % set chronux parameters if none provided
    if nargin < 5 || isempty(params)
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
        [p_part, t_part, f_part] = mtspecgramc(timetrial, [winsize winsize/4], params);
        
        % store
        power_tmp = 10*log10(permute(p_part, [4, 2, 3, 1]));
        % p_part: time x freqs x trials (x channels)
        powers(ch, :, :, :) = power_tmp; %#ok<SAGROW> % don't know how many freqs
        % powers: channels x frequencies x trials x time
        
        if ch==1
            freq = f_part';
            [~, I] = find_closest(freq, 300);
        end
        
    end
    
    fprintf('Calculating complete.\n')
    
    % place into struct, fix some vars
    data_out.trial = powers(:, 1:I, :, :);
    data_out.custom.filename = [data.custom.filename '_cmtsgrm'];
    data_out.custom.time = data.time{1}; % this is needed for some of the anova analysis
    
    % copy freqs
    data_out = rmfield(data_out, 'time');
    
    % fix data fields
    data_out.freq{1} = freq(1:I);
    data_out.freq_t = t_part;
    data_out.custom.params = params;
    
    % prepare data for save
    clear data;
    data = data_out; % used in save below
    save(data_out.custom.filename, 'data')
    clear data_out
    
end

fprintf('Exited function timecourse.\n')