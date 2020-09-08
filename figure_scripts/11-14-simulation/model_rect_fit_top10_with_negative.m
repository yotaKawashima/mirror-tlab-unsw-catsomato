%% Fit models to channels in S1 and S2
% First, find top 10% of channels that give the larger value of 
% sum(logSNR at harmonics and IMs).
% Then, fit models to the logSNR responses.
% Only Rectification
% Search coefficient from -1 to 1.

%% Set overall variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

%% Set script specific variables
data_dir = fullfile(data_path, 'included_datasets');

% find all cat names
cat_names = dirsinside(fullfile(data_path, 'included_datasets'));

% select datatype
datatype_logsnr = 'epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr';

% get date for saving the output files
dt = datestr(now, 'yyyymmdd');

%img_fmt = '-dpng';
img_fmt = '-depsc';

%% Harmonics and Intermodulations
f1 = 23;
f2 = 200;
f1_all = 23:23:250;
f1_harms = 46:23:250;
ims = sort([177:-23:0 200+23:23:250]);
harm_im = sort([f1_harms , ims]);     
foi_all = sort([f1_all, f2, ims]);            

%% Get logSNR from all channels (= bipolar channels x sessions) for each area
% Get only the max vibration condition from each session.
% Initialisation
data_struct = struct();

for area_id = 1:2
    area = ['S', num2str(area_id)]; % area
    
    % Initialisation 
    mean_logsnr_0_250_matrix = [];
    sd_logsnr_0_250_matrix = [];
    mean_logsnr_at_harms_ims_matrix = [];
    session_matrix = [];
    bipolar_chs_matrix =[];
    
    for session_id = 1:numel(cat_names)
        % find file names
        session_path = fullfile(data_dir, cat_names{session_id}, datatype_logsnr);
        loadname = dir(fullfile(session_path, ['*', area, '*.mat']));
        % load only the max vibration condition
        load(fullfile(session_path, loadname(end).name));
        
        % Get logSNR data
        % check if unipolar channels are included. If so, remove the
        % unipolar channels.
        if strcmp(data.label{1}(1:3), 'raw')
            % the first bipolar channel
            bipolarch_beginning = 1 + prod(data.custom.spatialconfig);
            % # of bipolar channels
            nChan = data.custom.nsignals - bipolarch_beginning + 1; 
        else
            bipolarch_beginning = 1;
            nChan = data.custom.nsignals;
        end
        
        % Get only bipolar channels data.
        % data dim = (bipolar channels x freqs x trials)
        logsnrs_ = data.trial(bipolarch_beginning:end, :, :);
        
        % Get logsnrs from 0Hz to 250Hz
        [~, fmax_ind] = find_closest(data.freq{1}, 250);
        logsnrs_0_250 = logsnrs_(:, 1:fmax_ind, :);
        
        % Mean across trials
        % data dim = bipolar channels x freqs
        logsnrs_0_250_mean_across_trials = squeeze(mean(logsnrs_0_250, 3));
        logsnrs_0_250_sd_across_trials = squeeze(std(logsnrs_0_250, 0, 3));
        
        % Get freq ids corresponding to harmonics and IMs
        [~, harm_im_inds] = find_closest(data.freq{1}, harm_im);
        
        % Get logSNR at harmonics and IMs
        % data dim = bipolar channels x (harmonics and ims)
        logsnrs_at_harm_im_mean_across_trials = ...
                logsnrs_0_250_mean_across_trials(:, harm_im_inds, :);
        
        % Store data across sessions
        % data dim = (bipolar channels x sessions) x (harmonics and ims)
        mean_logsnr_at_harms_ims_matrix = [mean_logsnr_at_harms_ims_matrix; ...
                                      logsnrs_at_harm_im_mean_across_trials];
        % data dim = (bipolar channels x sessions) x freqs
        mean_logsnr_0_250_matrix = [mean_logsnr_0_250_matrix;...
                                    logsnrs_0_250_mean_across_trials]; 
        sd_logsnr_0_250_matrix = [sd_logsnr_0_250_matrix;...
                                  logsnrs_0_250_sd_across_trials];
        session_matrix = [session_matrix;...
                          repmat(cat_names{session_id}, nChan, 1)];
        bipolar_chs_matrix = [bipolar_chs_matrix; (1:nChan)'];
    end % for session_id = 1:numel(cat_names)
    % Here, take sum of logsnrs across harms and ims. Based on the sum,
    % take top ??% channels.
    sum_logsnrs_matrix = sum(mean_logsnr_at_harms_ims_matrix, 2);
    
    % Store data with structure for each area
    data_struct(area_id).area = area;
    data_struct(area_id).freqs = data.freq{1}(1:fmax_ind);
    data_struct(area_id).session = session_matrix;
    data_struct(area_id).bipolar_ch = bipolar_chs_matrix;
    data_struct(area_id).sum_logsnrs_matrix = sum_logsnrs_matrix;
    data_struct(area_id).mean_logsnr_at_harms_ims_matrix = mean_logsnr_at_harms_ims_matrix;
    data_struct(area_id).mean_logsnr_0_250_matrix = mean_logsnr_0_250_matrix;
    data_struct(area_id).sd_logsnr_0_250_matrix = sd_logsnr_0_250_matrix;
    
end % for area_id = 1:2


%% Get thresholds for each area
% top 5% 
top10_threshold_S1 = prctile(data_struct(1).sum_logsnrs_matrix, 90, 1);
top10_threshold_S2 = prctile(data_struct(2).sum_logsnrs_matrix, 90, 1);

% Get top 5% channels for each area 
top10_channels_S1 = find(data_struct(1).sum_logsnrs_matrix > top10_threshold_S1);
top10_channels_S2 = find(data_struct(2).sum_logsnrs_matrix > top10_threshold_S2);

% Get log SNR of the 5% (mean across trials)
top10_mean_logsnrs_S1 = data_struct(1).mean_logsnr_0_250_matrix(top10_channels_S1, :);
top10_mean_logsnrs_S2 = data_struct(2).mean_logsnr_0_250_matrix(top10_channels_S2, :);

% Store data for each area
data_top10_struct = struct();
data_top10_struct(1).top_channels = top10_channels_S1;
data_top10_struct(2).top_channels = top10_channels_S2;
data_top10_struct(1).session = data_struct(1).session(top10_channels_S1, :);
data_top10_struct(2).session = data_struct(2).session(top10_channels_S2, :);
data_top10_struct(1).bipolar_ch = data_struct(1).bipolar_ch(top10_channels_S1, 1);
data_top10_struct(2).bipolar_ch = data_struct(2).bipolar_ch(top10_channels_S2, 1);
data_top10_struct(1).top_threshold = top10_threshold_S1;
data_top10_struct(2).top_threshold = top10_threshold_S2;
data_top10_struct(1).mean_logsnrs = top10_mean_logsnrs_S1;
data_top10_struct(2).mean_logsnrs = top10_mean_logsnrs_S2;

% Save data
save('Harmonics_IMs_top10_with_negative.mat', 'data_top10_struct');

%% Fitting some Model to actual data for comparison. 
% Try some models for the top 5% channels from S1 and S2.

%% Set variables
% Sampling rate [Hz]
sampling_rate = 10000; 

% time 0s to 2s.
time_s = 0; % Start [s]
time_e = 2; % End [s]
times = time_s:1/sampling_rate:time_e;

% Parameter for line search
tol = 1e-6;  
iteration = 20;
varmin = -1;
varmax = 1;


%% Fit model to logSNR for each area (S1 and S2)
each_area = struct();

for area_id = 1:2
    % load logSNR (mean across trials)
    logsnrs = data_top10_struct(area_id).mean_logsnrs;
    
    each_channel = struct();
    
    for bp_channel = 1:size(logsnrs, 1)
        logsnrs_now = logsnrs(bp_channel, :); % logSNR at this channel
        
        each_model_best = struct();
        
        for model_id = 1:4        
            switch model_id 
                case 1 % aRect(X) + bRect(Y)
                    model = @(a, b)model_plus(a, b);
                    % Line search
                    [L_min_var1, L_min_var2, L_error, L_errors, L_coordinate] = ...
                        fit_model_param2(model, varmin, varmax, logsnrs_now, sampling_rate, iteration, tol);
                    best_parameters = [L_min_var1, L_min_var2];
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2);

                case 2 % aRect(X) + bRect(Y) + cRect(X)Rect(Y)
                    model = @(a, b, c)model_without_RXY(a, b, c); 
                    % Line search
                    [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                        fit_model_param3(model, varmin, varmax, logsnrs_now, sampling_rate, iteration, tol);
                    best_parameters = [L_min_var1, L_min_var2, L_min_var3];
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2, L_min_var3);

                case 3 % aRect(X) + bRect(Y) + cRect(XY)
                    model = @(a, b, c)model_without_RXRY(a, b, c);
                    % Line search
                    [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                        fit_model_param3(model, varmin, varmax, logsnrs_now, sampling_rate, iteration, tol);
                    best_parameters = [L_min_var1, L_min_var2, L_min_var3];
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2, L_min_var3);

                case 4 % aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY) 
                    model = @(a, b, c, d)model_all(a, b, c, d); 
                    % Line search
                    [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                        fit_model_param4(model, varmin, varmax, logsnrs_now, sampling_rate, iteration, tol);
                    best_parameters = [L_min_var1, L_min_var2, ...
                                      L_min_var3, L_min_var4];
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
            end % switch model_id 
            
            % Compute logSNR for the best model
            [model_logsnrs, powers, freqs_] = compute_logsnrs_y(model_signal, sampling_rate);
            
            each_model_best(model_id).error = L_error;
            each_model_best(model_id).parameters = best_parameters;
            each_model_best(model_id).model_signal = model_signal;
            each_model_best(model_id).model_power = powers;
            each_model_best(model_id).model_logsnr = model_logsnrs;
            
        end % model_id = 1:4
        % Get the best among models
        [best_model_error, best_model_id] = min([each_model_best.error]);

        % Store data per 
        each_channel(bp_channel).logsnr = logsnrs_now;
        each_channel(bp_channel).each_model = each_model_best;
        each_channel(bp_channel).best_model_id = best_model_id;
        each_channel(bp_channel).best_model_error = best_model_error;        
        each_channel(bp_channel).best_model_signal = ...
                            each_model_best(best_model_id).model_signal;
        each_channel(bp_channel).best_model_logsnr = ...
                            each_model_best(best_model_id).model_logsnr;
    end % bp_channel
    
    each_area(area_id).area =  data_struct(area_id).area; % area
    each_area(area_id).freqs = data_struct(area_id).freqs; % frequencies
    each_area(area_id).each_channel = each_channel;
    
end % area_id = 1:2

% Save data
save('Results_Rect_top10.mat', 'each_area');

%% Functions
function [min_var1, min_var2, error, errors, coordinate] = fit_model_param2(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol)
% Find var1 and var 2 giving the local minimum by line serach. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  func = function to be accessed. Consider two-variable
%                  function (Model)
%           varmin = min of variable 
%           varmax = max of variable 
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
%           iteration = how many times you want to repeat linear search
%           tol = error tolerance for linear search
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%           coordinate = all var1 and var 2

%   Note:  Data should have the same sampling rate with func.
%   

    coordinate = nan(iteration, 2);
    errors = nan(iteration, 1);
    
    for i_iter = 1:iteration
        if mod(i_iter, 2) == 1 % Fit var1 first
            target = 1;
            if i_iter == 1
                varfixed = 0; % 
            else 
                varfixed = coordinate(i_iter - 1, 2);
            end % i_iter == 1
            
            [new_var, error]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [new_var, varfixed]; % new coordinate 
        elseif mod(i_iter, 2) == 0 % Then fit var2 
            target = 2;
            varfixed = coordinate(i_iter - 1, 1);
            [new_var, error]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
            coordinate_now = [varfixed, new_var]; % new coordinate            
        end % mod(i_iter, 2) == 1  
        
        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = error;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
end

function [min_var1, min_var2, min_var3, error, errors, coordinate] = fit_model_param3(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol)
% Find var1 and var 2 giving the local minimum by line serach. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  func = function to be accessed. Consider two-variable
%                  function (Model)
%           varmin = min of variable 
%           varmax = max of variable 
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
%           iteration = how many times you want to repeat linear search
%           tol = error tolerance for linear search
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
%           min_var3 = var3 giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%           coordinate = all var1, var2, and var3

%   Note:  Data should have the same sampling rate with func.
%   

    coordinate = nan(iteration, 3);
    errors = nan(iteration, 1);
    
    for i_iter = 1:iteration
        if mod(i_iter, 3) == 1 % Fit var1 first
            target = 1;
            if i_iter == 1
                varfixed = [0, 0]; 
            else 
                varfixed = coordinate(i_iter - 1, 2:3);
            end % i_iter == 1
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [new_var, varfixed(1), varfixed(2)]; % new coordinate 
        elseif mod(i_iter, 3) == 2 % Then fit var2 
            target = 2;
            if i_iter == 2
                varfixed = [coordinate(i_iter - 1, 1), 0]; 
            else 
                varfixed = [coordinate(i_iter - 1, 1), coordinate(i_iter - 1, 3)];
            end % i_iter == 1
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
            coordinate_now = [varfixed(1), new_var, varfixed(2)]; % new coordinate            
        elseif mod(i_iter, 3) == 0
            target = 3;
            varfixed = coordinate(i_iter - 1, 1:2);
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [varfixed(1), varfixed(2), new_var]; % new coordinate                         
        end % mod(i_iter, 2) == 1  
        
        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = error;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
    min_var3 = coordinate_now(3);
end

function [min_var1, min_var2, min_var3, min_var4, error, errors, coordinate] = fit_model_param4(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol)
% Find var1 and var 2 giving the local minimum by line serach. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  func = function to be accessed. Consider two-variable
%                  function (Model)
%           varmin = min of variable 
%           varmax = max of variable 
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
%           iteration = how many times you want to repeat linear search
%           tol = error tolerance for linear search
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
%           min_var3 = var3 giving the min difference
%           min_var4 = var4 giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%           coordinate = all var1, var2, and var3

%   Note:  Data should have the same sampling rate with func.
%   

    coordinate = nan(iteration, 4);
    errors = nan(iteration, 1);
    
    for i_iter = 1:iteration
        if mod(i_iter, 4) == 1 % Fit var1 first
            target = 1;
            if i_iter == 1
                varfixed = [0, 0, 0]; 
            else 
                varfixed = coordinate(i_iter - 1, 2:4);
            end % i_iter == 1
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [new_var, varfixed(1), varfixed(2), varfixed(3)]; % new coordinate 
        elseif mod(i_iter, 4) == 2 % Fit var2 
            target = 2;
            if i_iter == 2
                varfixed = [coordinate(i_iter - 1, 1), 0, 0]; 
            else 
                varfixed = [coordinate(i_iter - 1, 1), coordinate(i_iter - 1, 3:4)];
            end % i_iter == 1
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
            coordinate_now = [varfixed(1), new_var, varfixed(2), varfixed(3)]; % new coordinate            
        elseif mod(i_iter, 4) == 3 % Fit var3
            target = 3;
            if i_iter == 3
                varfixed = [coordinate(i_iter - 1, 1), coordinate(i_iter - 1, 2), 0];                
            else
                varfixed = [coordinate(i_iter - 1, 1:2), coordinate(i_iter - 1, 4)];
            end
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [varfixed(1), varfixed(2), new_var, varfixed(3)]; % new coordinate                         
        elseif mod(i_iter, 4) == 0
            target = 4;
            varfixed = coordinate(i_iter - 1, 1:3);
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [varfixed(1), varfixed(2), varfixed(3), new_var]; % new coordinate                                     
        end % mod(i_iter, 4) == 1  
        
        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = error;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
    min_var3 = coordinate_now(3);
    min_var4 = coordinate_now(4);
end

function [var, error]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
% Find var giving the closest value local minimum by golden section method.
%   Input:  model  = model to be accessed. Consider two-variable
%           varmin = min of variable 1
%           varmax = max of variable 1
%           varfixed = [fixed_var1]
%           datalogsnrs = data that we fit the model to.
%           sampling_rate = sampling rate for func
%           target = 1 (find var1) or 2 (find var2)
%           tol = error tolerance of variable
% 
%   Output: var = var giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%   Note:  Data should have the same sampling rate with func.
%   
    r = (sqrt(5) - 1)/2;

    % Set initial a,b,c,d and their func(x) values i.e. error(x)
    a = varmin;
    d = varmax;
    b = r*a + (1-r)*d;  % 1-r:r
    c = (1-r)*a + r*d;  % r:1-r
    
    if target == 1 % find var1 and fix var2
        signal_b = model(b, varfixed(1));
        signal_c = model(c, varfixed(1));
    elseif target == 2 % fix var1 and find var2
        signal_b = model(varfixed(1), b);
        signal_c = model(varfixed(1), c);
    end
    
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
     
    tol_ = d - a;

    % Update variabel until (d - a) < tol
    while(tol < tol_)
        if f_b < f_c % Eliminate x > c by replacing d with c.
            %a = a;
            d = c;
            c = b; 
            f_c = f_b; % func(c) can be copied from the pervious func(b).
            b = r*a + (1-r)*d;  % 1-r:r
            if target == 1 % update var1 and fix var2 
                signal_b = model(b, varfixed(1));  
            elseif target == 2 % fix var1 and update var2
                signal_b = model(varfixed(1), b);  
            end
            f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
            tol_ = d - a;
        else % Eliminate x < b by replacing a with b.
            a = b;
            %d = d;
            b = c;
            f_b = f_c; % func(d) can be copied from the pervious func(c).
            c = (1-r)*a + r*d;  % r:1-r
            if target == 1 % update var1 and fix var2 
                signal_c = model(c, varfixed(1));  
            elseif target == 2 % fix var1 and update var2
                signal_c = model(varfixed(1), c);  
            end           
            f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
            tol_ = d - a; 
        end % f_b < f_c
    end
    
    % Get the best var and its error
    if f_b < f_c
        var = b;
        error = f_b;
    else
        var = c;
        error = f_c;
    end % if f_b < f_c
    
end

function [var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
% Find var giving the closest value local minimum by golden section method.
%   Input:  model  = model to be accessed. Consider two-variable
%           varmin = min of variable 1
%           varmax = max of variable 1
%           varfixed = [fixed_var1, fixed_var2]
%           datalogsnrs = data that we fit the model to.
%           sampling_rate = sampling rate for func
%           target = 1 (find var1) or 2 (find var2) or 3 (find var3)
%           tol = error tolerance of variable
% 
%   Output: var = var giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%   Note:  Data should have the same sampling rate with func.
%   
    r = (sqrt(5) - 1)/2;

    % Set initial a,b,c,d and their func(x) values i.e. error(x)
    a = varmin;
    d = varmax;
    b = r*a + (1-r)*d;  % 1-r:r
    c = (1-r)*a + r*d;  % r:1-r
    
    if target == 1 % find var1 and fix var2 and var3
        signal_b = model(b, varfixed(1), varfixed(2));
        signal_c = model(c, varfixed(1), varfixed(2));
    elseif target == 2 % fix var1 and var3 and find var2
        signal_b = model(varfixed(1), b, varfixed(2));
        signal_c = model(varfixed(1), c, varfixed(2));
    elseif target == 3 % fix var1 and var2 and find var3
        signal_b = model(varfixed(1), varfixed(2), b);
        signal_c = model(varfixed(1), varfixed(2), c);       
    end
    
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
     
    tol_ = d - a;

    % Update variabel until (d - a) < tol
    while(tol < tol_)
        if f_b < f_c % Eliminate x > c by replacing d with c.
            %a = a;
            d = c;
            c = b; 
            f_c = f_b; % func(c) can be copied from the pervious func(b).
            b = r*a + (1-r)*d;  % 1-r:r
            if target == 1 % update var1 and fix var2 and var3
                signal_b = model(b, varfixed(1), varfixed(2));  
            elseif target == 2 % fix var1 and var3 and update var2
                signal_b = model(varfixed(1), b, varfixed(2));  
            elseif target == 3 % fix var1 and var2 and update var3
                signal_b = model(varfixed(1), varfixed(2), b);
            end
            f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
            tol_ = d - a;
        else % Eliminate x < b by replacing a with b.
            a = b;
            %d = d;
            b = c;
            f_b = f_c; % func(d) can be copied from the pervious func(c).
            c = (1-r)*a + r*d;  % r:1-r
            if target == 1 % update var1 and fix var2 and var3
                signal_c = model(c, varfixed(1), varfixed(2));  
            elseif target == 2 % fix var1 and var3 and update var2
                signal_c = model(varfixed(1), c, varfixed(2));  
            elseif target == 3 % fix var1 and var2 and update var3
                signal_c = model(varfixed(1), varfixed(2), c);                
            end           
            f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
            tol_ = d - a; 
        end % f_b < f_c
    end
    
    % Get the best var and its error
    if f_b < f_c
        var = b;
        error = f_b;
    else
        var = c;
        error = f_c;
    end % if f_b < f_c
    
end

function [var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
% Find var giving the closest value local minimum by golden section method.
%   Input:  model  = model to be accessed. Consider two-variable
%           varmin = min of variable 1
%           varmax = max of variable 1
%           varfixed = [fixed_var1, fixed_var2, fixed_var3]
%           datalogsnrs = data that we fit the model to.
%           sampling_rate = sampling rate for func
%           target = 1 (find var1) or 2 (find var2) or 3 (find var3) or
%                    4 (find var4)
%           tol = error tolerance of variable
% 
%   Output: var = var giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%   Note:  Data should have the same sampling rate with func.
%   
    r = (sqrt(5) - 1)/2;

    % Set initial a,b,c,d and their func(x) values i.e. error(x)
    a = varmin;
    d = varmax;
    b = r*a + (1-r)*d;  % 1-r:r
    c = (1-r)*a + r*d;  % r:1-r
    
    if target == 1 % find var1 and fix var2, var3, and var4
        signal_b = model(b, varfixed(1), varfixed(2), varfixed(3));
        signal_c = model(c, varfixed(1), varfixed(2), varfixed(3));
    elseif target == 2 % fix var1, var3, and var4 and find var2
        signal_b = model(varfixed(1), b, varfixed(2), varfixed(3));
        signal_c = model(varfixed(1), c, varfixed(2), varfixed(3));
    elseif target == 3 % fix var1, var2, and var4 and find var3
        signal_b = model(varfixed(1), varfixed(2), b, varfixed(3));
        signal_c = model(varfixed(1), varfixed(2), c, varfixed(3));       
    elseif target == 4 % fix var1, var2, and var3 and find var4
        signal_b = model(varfixed(1), varfixed(2), varfixed(3), b);
        signal_c = model(varfixed(1), varfixed(2), varfixed(3), c);        
    end
    
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
     
    tol_ = d - a;

    % Update variabel until (d - a) < tol
    while(tol < tol_)
        if f_b < f_c % Eliminate x > c by replacing d with c.
            %a = a;
            d = c;
            c = b; 
            f_c = f_b; % func(c) can be copied from the pervious func(b).
            b = r*a + (1-r)*d;  % 1-r:r
            if target == 1 % update var1 and fix var2, var3, and var4
                signal_b = model(b, varfixed(1), varfixed(2), varfixed(3));  
            elseif target == 2 % fix var1, var3, and var4 and update var2
                signal_b = model(varfixed(1), b, varfixed(2), varfixed(3));  
            elseif target == 3 % fix var1, var2, and var4 and update var3
                signal_b = model(varfixed(1), varfixed(2), b, varfixed(3));
            elseif target == 4 % fix var1, var2, and var3 and update var4
                signal_b = model(varfixed(1), varfixed(2), varfixed(3), b);
            end
            f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
            tol_ = d - a;
        else % Eliminate x < b by replacing a with b.
            a = b;
            %d = d;
            b = c;
            f_b = f_c; % func(d) can be copied from the pervious func(c).
            c = (1-r)*a + r*d;  % r:1-r
            if target == 1 % update var1 and fix var2, var3, and var4
                signal_c = model(c, varfixed(1), varfixed(2), varfixed(3));  
            elseif target == 2 % fix var1, var3, and var4 and update var2
                signal_c = model(varfixed(1), c, varfixed(2), varfixed(3));  
            elseif target == 3 % fix var1, var2, and var4 and update var3
                signal_c = model(varfixed(1), varfixed(2), c, varfixed(3));                
            elseif target == 4 % fix var1, var2, and var3 and update var4
                signal_c = model(varfixed(1), varfixed(2), varfixed(3), c);
            end           
            f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
            tol_ = d - a; 
        end % f_b < f_c
    end
    
    % Get the best var and its error
    if f_b < f_c
        var = b;
        error = f_b;
    else
        var = c;
        error = f_c;
    end % if f_b < f_c
    
end

function error = compute_error(model_signal, datalogsnrs, sampling_rate)
% Compute error
%   Input: sampling_rate
%   Output: error

    % Frequencies of interest (Just for visualisatoin)
    foi_f1_and_harm = 23:23:250;
    foi_f2 = 200;
    foi_inter = sort([177:-23:0 200+23:23:250]);
    foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);            

    % Compute logSNR
    [modellogsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_signal, sampling_rate);
    
    [~, foi_inds] = find_closest(freqs, foi);
    error = sum(abs(modellogsnrs(1, foi_inds) - datalogsnrs(1, foi_inds)));
    
    % Check how log power looks like
    %{
    figure();clf
    [~, maxf_ind] = find_closest(freqs, 250);
    plot(freqs(1:maxf_ind), powers_tmp(1:maxf_ind))
    close gcf;    
    %}
    
end

function [logsnrs, powers, freqs] = compute_logsnrs_y(data, sampling_rate)
% Compute power spectrum and return logSNR.
% The same logsnr computation as Renee's snrtosurrounds.m, but this code is
% more easy to follow what we are doing.
%
%   The neighbourhood region is illustrated in the following diagram.
%   F=frequency of interest
%   A1=F - outerdiff       F - 3 Hz
%   B1=F - innerdiff       F - 1 Hz
%   C1=F - halfbandwidth   F - 0.5Hz
%   C2=F + halfbandwidth   F + 0.5Hz
%   B2=F + innerdiff       F + 1 Hz
%   A2=F + outerdiff       F + 3 Hz
%   
%   L = -1/(# of samples in low region)
%   H = +1/(# of samples in high region)
%   
%      LLLLLLLLLLL000HHHHH000LLLLLLLLLLL    <- Mask 
%   ... ......... ... ... ... ... ......... ..  <- data points 
%      |         |   |   |   |   |         |
%      A1        B1 C1   F   C2  B2        A2   <-freq
%      |         |   |   |   |   |         |
%     -3        -1  -0.5 0 +0.5  +1        +3
%                     ___ ___              
%    __           ___|       |___           __ 
%      |_________|               |_________|
%          low         high          low
%       
%   (This is just visualisation. # of points is not the same
%   as real data.

% Written by Yota. 
    % Parameters given to chronux mtspectrumc.m.
    params.Fs = sampling_rate;
    params.tapers = [1 1]; % hbw = (k+1)/2T  k=# of tapers, T=length of time
    params.pad = 1;   
    % Here, k = 1, T = 2s. Then, hbw = (1+1)/2*2 = 0.5Hz 
    halfbandwidth = 0.5; 
    
    % Compute power spectrum
    [p_part, f_part] = mtspectrumc(data, params);
    powers = 10*log10(p_part');
    freqs = f_part';

    % Parameter for computing logSNR
    innerdiff = 1;
    outerdiff = 3;
    half_width_low = outerdiff - innerdiff;
    maxf = 250;

    % create the mask
    % first, find the number of indicies I need to be 1 Hz
    % assume they're all evenly spaced. 
    spacing = freqs(2) - freqs(1);
    nSamp = round(1/spacing);    % # of indicies within 1Hz width. 
    nSamp_hbw = floor(nSamp*halfbandwidth); % # of indicies in hbw
    half_nSamp_low = half_width_low*nSamp; % # of indicies within low range
    
    % mask covering from A1 to A2.
    mLen = 2*(2*nSamp_hbw+half_nSamp_low)+1; % This could be wrong
    mask = zeros(1, mLen);

    % neighbour region outside C1-F-C2
    % -1/(# of Samples in neighbor region) 
    nSamp_low = 2*half_nSamp_low;      % # of indices in low
    lomask = 1/nSamp_low;                    % -1/(# of samples in low)
    mask(1:half_nSamp_low) = -lomask;         % from A1 to B1
    mask(end-half_nSamp_low+1:end) = -lomask; % from B2 to A2

    % central region C1-F-C2 (i.e. tagged freq)
    % # of samples in high 
    nSamp_high = nSamp_hbw*2 + 1; % Here, 1 is for data point of F
    himask = 1/nSamp_high;        %  +1/(# of samples in high) 
    m = ceil(mLen/2);             % Get the center of data point.
    mask(m-nSamp_hbw:m+nSamp_hbw) = himask;
    
    % reduce to only frequencies of interest
    [~, maxf_ind] = find_closest(freqs, maxf);
    freqs_interest = freqs(1:maxf_ind);
    
    % first frequency index at hbw*nS 
    % Some data points in the beginning and the end should be removed
    % because neighboring does not exist for them.
    % We can do zero-padding (to signal not to mask) to solve this. 
    % first frequency index at 3*nSample. 
    freqs = freqs_interest(nSamp_hbw*nSamp:(end-nSamp_hbw*nSamp+1));
    nReps = length(freqs);
    
    logsnrs = zeros(1, nReps);

    for index = 1:nReps
        % zero pad mask
        maskm = [zeros(1, index-1) mask zeros(1,maxf_ind-mLen-index+1)];
        logsnrs(1, index) = dot(maskm, powers(1:maxf_ind));
    end
    
    % Extract power within a range for 0Hz to 250Hz
    powers_interest = powers(1:maxf_ind);
    powers = powers_interest(nSamp_hbw*nSamp:(end-nSamp_hbw*nSamp+1));
    
end

%% Model 
function model = model_plus(a, b)
% aRect(X) + bRect(Y)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = (a*Rect_comp1) + (b*Rect_comp2);
end

function model = model_without_RXY(a, b, c)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
        
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_comp1 .* Rect_comp2;
end

function model = model_without_RXRY(a, b, c)
% aRect(X) + bRect(Y) + cRect(XY)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Rect_mult_sines = mult_sines;
    Rect_mult_sines(mult_sines <=0 ) = 0;
    
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_mult_sines;

end

function model = model_all(a, b, c, d)
% aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY)
    % frequency
    f1 = 23;
    f2 = 200;

    % Sampling rate [Hz]
    sampling_rate = 10000; 

    % time 0s to 2s.
    time_s = 0; % Start [s]
    time_e = 2; % End [s]
    times = time_s:1/sampling_rate:time_e;
    
    % Model
    % Components of models.
    % f1[Hz] sine wave.
    comp1 = sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    
    % multiply sine waves
    mult_sines = comp1 .* comp2;
    Rect_mult_sines = mult_sines;
    Rect_mult_sines(mult_sines <=0 ) = 0;
    
    % model 
    model = a*Rect_comp1 + b*Rect_comp2 + c*Rect_comp1 .* Rect_comp2 ...
            + d*Rect_mult_sines;
end
