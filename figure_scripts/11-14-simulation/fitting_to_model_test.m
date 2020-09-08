%% Fitting some Model to actual data for comparison. (Test)
% Here, we fit some model to actual data. Then, based on the we check
% spatial distribution of the fitted models for each area (S1 and S2).

%% Set variables
[fpath, fname, fext] = fileparts(mfilename('fullpath'));
run(fullfile(fpath, '../path_setup.m'));

catname = 'C20110808_R03';
area = 'S1';
% If you included unipolar channels in your stored data. You need to be
% careful that bipolar channel id starts at 1011 and at 65 for S1 and S2.
bp_chan_1 = 256; % bipolar channel
bp_chan_2 = 152; % bipolar channel

% Frequencies of interest (Just for visualisatoin)
foi_f1_and_harm = 23:23:250;
foi_f2 = 200;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi = sort([foi_f1_and_harm, foi_f2, foi_inter]);

% frequency
f1 = 23;
f2 = 200;

% Sampling rate [Hz]
sampling_rate = 10000; 

% time 0s to 2s.
time_s = 0; % Start [s]
time_e = 2; % End [s]
times = time_s:1/sampling_rate:time_e;

% parameter for model
a_range = [0:0.01:0.1];
b_range = [0:0.1:5];

%% Model
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

%model = @(a,b)a*Rect_comp1 + b*Rect_comp2;
%model = @(a,b)(a*Rect_comp1) .* (b*Rect_comp2);
model = @(a, b)model_plus(a,b); % a and b are inslide the Rect()
%model = @(a, b)model_mult(a,b); % a and b are inslide the Rect()

%%  Set data path to logSNR data that we fit a model to.
data_dir_logsnr = [data_path 'included_datasets/' catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'];
data_dir_bp = [data_path 'included_datasets/' catname '/epoched_rsampsl_biprref_evkresp_cmtspwr_snrsurr'];

%% Get intermodulation prominent bipolar channels. 
% Compure mean and std across trials for this channel data. 
[mean_logsnr_1, std_logsnr_1, freq_data_1] = mean_std_data(data_dir_logsnr, catname, bp_chan_1);

%% Rough grid search
[G_min_var1_1, G_min_var2_1, G_error_1, G_errors_1] = fit_model_v1(model, a_range, b_range, mean_logsnr_1, sampling_rate);

% Create grid for plot
coordinate_1 = nan(length(a_range)* length(b_range), 2);
increment_ = 1;
for a_var = a_range
    for b_var = b_range
        coordinate_1(increment_, :) =  [a_var, b_var];
        increment_ = increment_ + 1;
    end
end 
% Reshape errors dim = a_range x b_range -> (a_range*b_range) x 1
errors_1_resh = reshape(G_errors_1', [], 1);
figure();clf
scatter(coordinate_1(:, 1), coordinate_1(:, 2), [] ,errors_1_resh, 'filled');
cbar = colorbar();
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
cbar.Label.String = 'Error at foi';
set(gca, 'fontsize', 16);
figure();clf
plot3(coordinate_1(:, 1), coordinate_1(:, 2), errors_1_resh,'.','MarkerSize',15);
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
zlabel('Error at foi');
set(gca, 'fontsize', 16);

% Check results 
model_best_1 = model(G_min_var1_1, G_min_var2_1);
% Compute logSNR
[model_logsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_best_1, sampling_rate);
    
% Plot logSNR
figure();clf
hold on
plot(freqs, model_logsnrs);
plot(freqs, mean_logsnr_1);

% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(2)-ylims(1))*0.1+ylims(1)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 2, 'LineStyle', '--')
hold off
set(gca, 'XLim', [0,250]);
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
legend({'model', 'real data'});
%title_text = ['a*Rect(X) + b*Rect(Y)', ', a=', num2str(G_min_var1_1), ', b=', num2str(G_min_var2_1), ...
%               ', error=', num2str(G_error_1)];
title_text = ['a*Rect(X) * b*Rect(Y)', ', a=', num2str(G_min_var1_1), ', b=', num2str(G_min_var2_1), ...
               ', error=', num2str(G_error_1)];
title(title_text);

%% Line search
tol = 1e-6;  
iteration = 20;
varmin = 0;
varmax = 0.05;
[L_min_var1_1, L_min_var2_1, L_error_1, L_errors_1, L_coordinate_1] = fit_model_v3(model, varmin, varmax, mean_logsnr_1, sampling_rate, iteration, tol);

% Check 
figure();clf
scatter(L_coordinate_1(:, 1), L_coordinate_1(:, 2), [] ,L_errors_1, 'filled');
cbar = colorbar();
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
cbar.Label.String = 'Error at foi';
set(gca, 'fontsize', 16);
figure();clf
plot3(L_coordinate_1(:, 1), L_coordinate_1(:, 2), L_errors_1,'.','MarkerSize',15);
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
zlabel('Error at foi');
set(gca, 'fontsize', 16);

% Check results 
model_best_1 = model(L_min_var1_1, L_min_var2_1);
% Compute logSNR
[model_logsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_best_1, sampling_rate);
    
% Plot logSNR
figure();clf
hold on
plot(freqs, model_logsnrs);
plot(freqs, mean_logsnr_1);

% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(2)-ylims(1))*0.1+ylims(1)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 2, 'LineStyle', '--')
hold off
set(gca, 'XLim', [0,250]);
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
legend({'model', 'real data'});

%title_text = ['a*Rect(X) + b*Rect(Y)', ', a=', num2str(L_min_var1_1), ',b=', num2str(L_min_var2_1),...
%               ', error=', num2str(L_error_1)];
title_text = ['a*Rect(X) * b*Rect(Y)', ', a=', num2str(L_min_var1_1), ', b=', num2str(L_min_var2_1),...
              ', error=', num2str(L_error_1)];
title(title_text);

%% Get intermodulation non-prominent bipolar channels.

% Compure mean and std across trials for this channel data. 
[mean_logsnr_2, std_logsnr_2, freq_data_2] = mean_std_data(data_dir_logsnr, catname, bp_chan_2);

%% Rough Grid search 
[G_min_var1_2, G_min_var2_2, G_error_2, G_errors_2] = fit_model_v1(model, a_range, b_range, mean_logsnr_2, sampling_rate);

% Create grid for plot
coordinate_2 = nan(length(a_range)* length(b_range), 2);
increment_ = 1;
for a_var = a_range
    for b_var = b_range
        coordinate_2(increment_, :) =  [a_var, b_var];
        increment_ = increment_ + 1;
    end
end 

% Reshape errors dim = a_range x b_range -> (a_range*b_range) x 1
errors_2_resh = reshape(G_errors_2', [], 1);
figure();clf
scatter(coordinate_2(:, 1), coordinate_2(:, 2), [] ,errors_2_resh, 'filled');
cbar = colorbar();
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
cbar.Label.String = 'Error at foi';
set(gca, 'fontsize', 16);
figure();clf
plot3(coordinate_2(:, 1), coordinate_2(:, 2), errors_2_resh,'.','MarkerSize',15);
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
zlabel('Error at foi');
set(gca, 'fontsize', 16);

% Check results 
model_best_2 = model(G_min_var1_2, G_min_var2_2);

% Compute logSNR
[model_logsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_best_2, sampling_rate);


figure();clf
hold on
plot(freqs, model_logsnrs);
plot(freqs, mean_logsnr_2);

% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(2)-ylims(1))*0.1+ylims(1)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 2, 'LineStyle', '--')
hold off
set(gca, 'XLim', [0,250]);
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
legend({'model', 'real data'});
%title_text = ['a*Rect(X) + b*Rect(Y)', ', a=', num2str(G_min_var1_2), ', b=', num2str(G_min_var2_2),...
%               ', error=', num2str(G_error_2)];
title_text = ['a*Rect(X) * b*Rect(Y)', ', a=', num2str(G_min_var1_2), ', b=', num2str(G_min_var2_2),...
              ', error=', num2str(G_error_2)];
title(title_text);

%% Line search
[L_min_var1_2, L_min_var2_2, L_error_2, L_errors_2, L_coordinate_2] = fit_model_v3(model, varmin, varmax, mean_logsnr_2, sampling_rate, iteration, tol);

% Check 
figure();clf
scatter(L_coordinate_2(:, 1), L_coordinate_2(:, 2), [] ,L_errors_2, 'filled');
cbar = colorbar();
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
cbar.Label.String = 'Error at foi';
set(gca, 'fontsize', 16);
figure();clf
plot3(L_coordinate_2(:, 1), L_coordinate_2(:, 2), L_errors_2,'.','MarkerSize',15);
xlabel('a = amp of 23Hz');
ylabel('b = amp of 200Hz');
zlabel('Error at foi');
set(gca, 'fontsize', 16);

% Check results 
model_best_2 = model(L_min_var1_2, L_min_var2_2);
% Compute logSNR
[model_logsnrs, powers_tmp, freqs] = compute_logsnrs_y(model_best_2, sampling_rate);
    
% Plot logSNR
figure();clf
hold on
plot(freqs, model_logsnrs);
plot(freqs, mean_logsnr_2);

% plot vertical lines at frequencies of interest
% grab axis variables
ylims = get(gca, 'YLim');
v_ylims = [ylims(1), (ylims(2)-ylims(1))*0.1+ylims(1)];
v1 = numel(foi_f1_and_harm);
v2 = numel(foi_inter);
v = zeros(1, numel(foi));
% plot vertical lines
v(1:v1) = plot([foi_f1_and_harm', foi_f1_and_harm'], v_ylims, 'Color', [255, 47, 64]/256); % 23 + harmonics
v(v1+1) = plot([foi_f2, foi_f2], v_ylims, 'Color', [0, 205, 0]/256);%[83, 215, 81]/256); % 200
v(v1+2:v1+v2+1) = plot([sort(foi_inter)', sort(foi_inter)'], v_ylims, 'Color', [0, 82, 255]/256); % intermodulation
% format vertical lines
set(v, 'LineWidth', 2, 'LineStyle', '--')
hold off
set(gca, 'XLim', [0,250]);
xlabel('frequency [Hz]');
ylabel('logSNR [dB]');
legend({'model', 'real data'});
%title_text = ['a*Rect(X) + b*Rect(Y)', ', a=', num2str(L_min_var1_2), ', b=', num2str(L_min_var2_2),...
%               ', error=', num2str(L_error_1)];
title_text = ['a*Rect(X) * b*Rect(Y)', ', a=', num2str(L_min_var1_2), ', b=', num2str(L_min_var2_2),...
              ', error=', num2str(L_error_2)];

title(title_text);

%% Functions
% Compute mean across trials
function [mean_data, std_data, freq_data] = mean_std_data(data_dir, catname, chan)
% Load S1 data only max condition. Then compute mean and std across trials.
%   Input:  data_dir = directory name
%           fig1_catname = cat name 
%           chan = channels we plot
%   Output: mean_data = mean across trials
%           std_data = std across trials
%           freq_data = frequencies (0-250Hz)

    % want the max stimulation condition
    fnames = dir(fullfile(data_dir, [catname, '*S1*.mat'])); 
    nConds = length(fnames);
    % load the last one (max stimulation condition)
    load(fullfile(data_dir, fnames(nConds).name))

    % Extract data from the channels (dim = frequencies x trials)
    trial_data = squeeze(data.trial(chan,:,:)); 

    % Get index of 250Hz
    [~,fend] = find_closest(data.freq{1}, 250);

    % Take mean of log SNR across trials.
    mean_data = mean(trial_data(1:fend, :), 2);
    std_data = std(trial_data(1:fend, :), 0, 2); % unbiased
    % Correct dim
    mean_data = mean_data';
    std_data = std_data';
    % Frequency data 
    freq_data = data.freq{1}(1:fend);
    freq_data = freq_data';
end


function [min_var1, min_var2, error, errors] = fit_model_v1(model, var1range, var2range, datalogsnrs, sampling_rate)
% Find var1 and var 2 giving the local minimum by grid serach. We minimise
% error (difference of logSNR between model and data across foi).
%   Input:  model = model to be accessed. Consider two-variable
%           var1range = domain of function (variable 1)
%           var2range = domain of function (variable 2)
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
%           error = The smallest difference between function and data.
%           errors = All errors
%   Note:  Data should have the same sampling rate with func.
%   
    errors = nan(length(var1range), length(var2range));
    
    for var1_id = 1:length(var1range)
        var1 = var1range(var1_id);
        for var2_id = 1:length(var2range)
            var2 = var2range(var2_id);
            %disp(var1)
            %disp(var2)
            %disp('-------')
            % Compute error 
            if var1 == 0 && var2 == 0 % If 0 amp signal
                error_tmp = Inf;
            else 
                error_tmp = compute_error(model(var1, var2), datalogsnrs, sampling_rate);
            end
            
            % Store error 
            errors(var1_id, var2_id) = error_tmp;
            
            % Keep the first parameters and error
            % note that the parameters (0,0) results in logSNR of NaN 
            % because the wave form is line with 0 value. 
            if var1 == 0 && var2 == 0
                min_var1 = var1;
                min_var2 = var2;
                error = inf; % 
            end % var1 == var1range(1) && var2 == var2range(2)
            
            % Store parameters if error_tmp is smaller than the stored one.
            if error_tmp < error 
                min_var1 = var1;
                min_var2 = var2;
                error = error_tmp;         
            end % error_tmp < error    
            
        end % for var2 = var2range
    end % for var1 = var1range
end

%{
function [min_var1, min_var2, error, errors] = fit_model_v2(model, var1max, var2max, data, sampling_rate)
% Find var1 and var 2 giving the local minimum by Nelder-Mead method. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  model = model to be accessed. Consider two-variable
%           var1max = max of variable 1
%           var2max = max of variable 2
%           data = data that we fit the function to.
%           sampling_rate = sampling rate for func
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
%           error = The smallest difference between function and data.
%           errors = all errors
% 
%   Note:  Data should have the same sampling rate with func.
% 

    % Define the initial simplex as tabel
    coordinate = [0, 0; var1max, 0; 0, var2max];
    error1 = inf;
    error2 = compute_error(model(var1max, 0), data, sampling_rate);
    error3 = compute_error(model(0, var2max), data, sampling_rate); 
    error = [error1; error2; error3];
    simplex = table(coordinate, error);
    
    errors = nans(1000,1);

    % Update simplex based on error 
    % Iteration stops at 100 (No other termination criteria)    
    for i=1:1000
        % Sort vertices of the simplex based on error 
        % (the first is the smallest error).
        simplex = sortrows(simplex, 'error');  
        
        errors(i) = simplex.error(1);
        
        % Compute the centroid
        centroid = (simplex.coordinate(1, :) + simplex.coordinate(2, :))/2;

        % Reflection point
        refpoint = 2*centroid - simplex.coordinate(3, :);
        f_refpoint = compute_error(model(refpoint(1), refpoint(2)), ...
                                   data, sampling_rate);
                               
        % If f(p1) <= f(ref) <= f(pm)
        if (simplex.error(1) <= f_refpoint)&&(f_refpoint <= simplex.error(2))
            if (refpoint(1) >= 0)&&(refpoint(2) >= 0)
                simplex.coordinate(3, :) = refpoint;
                simplex.error(3) = f_refpoint;
                continue
            else 
                % Reduction
                simplex.coordinate(2, :) = 1.5*simplex.coordinate(2, :) ...
                                          - 0.5*simplex.coordinate(1, :);
                simplex.error(2) = compute_error(model(simplex.coordinate(2, 1), ...
                                            simplex.coordinate(2, 2)), ...
                                            data, sampling_rate);
                simplex.coordinate(3, :) = 1.5*simplex.coordinate(3, :) ...
                                          - 0.5*simplex.coordinate(1, :);
                simplex.error(3) = compute_error(model(simplex.coordinate(3, 1), ...
                                            simplex.coordinate(3, 2)), ...
                                            data, sampling_rate);
                continue
            end
          
        % If f(ref) < f(p1)    
        elseif f_refpoint < simplex.error(1)
            % Expansion point
            exppoint = 3*centroid - 2*simplex.coordinate(3, :);
            f_exppoint = compute_error(model(exppoint(1), exppoint(2)), ...
                                       data, sampling_rate);
            % if f(exp) < f(ref) < f(p1)
            if (f_exppoint < f_refpoint)&&(f_refpoint < simplex.error(1))
                if (exppoint(1) >= 0)&&(exppoint(2) >= 0)
                    simplex.coordinate(3, :) = exppoint;
                    simplex.error(3) = f_exppoint;
                else
                   % Reduction
                   simplex.coordinate(2, :) = 1.5*simplex.coordinate(2, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                   simplex.error(2) = compute_error(model(simplex.coordinate(2, 1), ...
                                                simplex.coordinate(2, 2)), ...
                                                data, sampling_rate);
                   simplex.coordinate(3, :) = 1.5*simplex.coordinate(3, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                   simplex.error(3) = compute_error(model(simplex.coordinate(3, 1), ...
                                                simplex.coordinate(3, 2)), ...
                                                data, sampling_rate);                    
                end
            % if f(ref) <= f(exp)
            elseif f_refpoint <= f_exppoint
                if (refpoint(1) >= 0)&&(refpoint(2) >= 0)
                    simplex.coordinate(3, :) = refpoint;
                    simplex.error(3) = f_refpoint;
                    continue
                else 
                    % Reduction
                    simplex.coordinate(2, :) = 1.5*simplex.coordinate(2, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                    simplex.error(2) = compute_error(model(simplex.coordinate(2, 1), ...
                                                simplex.coordinate(2, 2)), ...
                                                data, sampling_rate);
                    simplex.coordinate(3, :) = 1.5*simplex.coordinate(3, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                    simplex.error(3) = compute_error(model(simplex.coordinate(3, 1), ...
                                                simplex.coordinate(3, 2)), ...
                                                data, sampling_rate);
                    continue
                end
            end %(f_exppoint < f_refpoint)&&(f_refpoint < simplex{2,1})
        
        % If f(m) <= f(ref)
        elseif simplex.error(2) <= f_refpoint
            % Contraction point
            conpoint = 1/2*centroid + 1/2*simplex.coordinate(3, :);
            f_conpoint = compute_error(model(conpoint(1), conpoint(2)), ...
                                       data, sampling_rate);
            % if f(cont) <= f(pm+1)
            if f_conpoint <= simplex.error(3)
                if (conpoint(1) >= 0)&&(conpoint(2) >= 0)
                    simplex.coordinate(3, :) = conpoint;
                    simplex.error(3) = f_conpoint;
                    continue
                else 
                    % Reduction
                    simplex.coordinate(2, :) = 1.5*simplex.coordinate(2, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                    simplex.error(2) = compute_error(model(simplex.coordinate(2, 1), ...
                                                simplex.coordinate(2, 2)), ...
                                                data, sampling_rate);
                    simplex.coordinate(3, :) = 1.5*simplex.coordinate(3, :) ...
                                              - 0.5*simplex.coordinate(1, :);
                    simplex.error(3) = compute_error(model(simplex.coordinate(3, 1), ...
                                                simplex.coordinate(3, 2)), ...
                                                data, sampling_rate);
                    continue
                end
            % if f(pm+1) < f(cont)
            elseif simplex.error(3) < f_conpoint
               % Reduction
               simplex.coordinate(2, :) = 1.5*simplex.coordinate(2, :) ...
                                          - 0.5*simplex.coordinate(1, :);
               simplex.error(2) = compute_error(model(simplex.coordinate(2, 1), ...
                                            simplex.coordinate(2, 2)), ...
                                            data, sampling_rate);
               simplex.coordinate(3, :) = 1.5*simplex.coordinate(3, :) ...
                                          - 0.5*simplex.coordinate(1, :);
               simplex.error(3) = compute_error(model(simplex.coordinate(3, 1), ...
                                            simplex.coordinate(3, 2)), ...
                                            data, sampling_rate);
            end % if f(cont) <= f(pm+1)
        end % if simplex{2,1} < f_refpoint && simplex{2,2} < f_refpoint    
    end % i=1:100 
    
    % Get the best vertix from the last simplex
    simplex = sortrows(simplex, 'error');  
    min_var1 = simplex.coordinate(1, 1); 
    min_var2 = simplex.coordinate(1, 2);
    error = simplex.error(1);
end
%}

function [min_var1, min_var2, error, errors, coordinate] = fit_model_v3(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol)
% Find var1 and var 2 giving the local minimum by linear serach. We
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
            
            [new_var, error]= golden_section(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
            coordinate_now = [new_var, varfixed]; % new coordinate 
        elseif mod(i_iter, 2) == 0 % Then fit var2 
            target = 2;
            varfixed = coordinate(i_iter - 1, 1);
            [new_var, error]= golden_section(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
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

function [var, error]= golden_section(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
% Find var giving the closest value local minimum by golden section method.
%   Input:  model  = model to be accessed. Consider two-variable
%           varmin = min of variable 1
%           varmax = max of variable 1
%           varfixed = value of fixed variable
%           datalogsnrs = data that we fit the model to.
%           sampling_rate = sampling rate for func
%           target = 1 (find var1, fix var2) or 2 (find var2, fix var1)
%           tol = error tolerance of variable
% 
%   Output: min_var1 = var1 giving the min difference
%           min_var2 = var2 giving the min difference
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
        signal_b = model(b, varfixed);
        signal_c = model(c, varfixed);
    elseif target == 2 % fix var1 and find var2
        signal_b = model(varfixed, b);
        signal_c = model(varfixed, c);
    end
    
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
     
    tol_ = d - a;

    % Update variabel until (d - a) < tol
    while(tol < tol_)
        figure();clf
        scatter([a, b, c, d], [1, 1, 1, 1]);
        close all;
        if f_b < f_c % Eliminate x > c by replacing d with c.
            %a = a;
            d = c;
            c = b; 
            f_c = f_b; % func(c) can be copied from the pervious func(b).
            b = r*a + (1-r)*d;  % 1-r:r
            if target == 1 % update var1 and fix var2
                signal_b = model(b, varfixed);
            elseif target == 2 % fix var1 and update var2
                signal_b = model(varfixed, b);
            end
            f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
            tol_ = d - a;
        else % Eliminate x < b by replacing a with b.
            a = b;
            %d = d;
            b = c;
            f_b = f_c; % func(d) can be copied from the pervious func(c).
            c = (1-r)*a + r*d;  % r:1-r
            if target == 1 % find var1 and fix var2
                signal_c = model(c, varfixed);
            elseif target == 2 % fix var1 and find var2
                signal_c = model(varfixed, c);
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

end

function model = model_plus(a, b)
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
    comp1 = a*sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = b*sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = Rect_comp1 + Rect_comp2;
end

function model = model_mult(a, b)
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
    comp1 = a*sin(2*pi*f1*times);

    % f2[Hz] sine wave.
    comp2 = b*sin(2*pi*f2*times);

    % Rectification 
    Rect_comp1 = comp1;
    Rect_comp1(comp1 <= 0) = 0;
    Rect_comp2 = comp2;
    Rect_comp2(comp2 <= 0) = 0;
    model = Rect_comp1 .* Rect_comp2;
end
