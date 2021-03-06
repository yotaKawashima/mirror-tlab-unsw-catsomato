%% Plot exemplar responses (Test)
% First, find top 10% of channels that give the larger value of 
% sum(logSNR at harmonics and IMs).
% Then, fit models to the logSNR responses.
% Only Rectification
% Search coefficient from -1 to 1.
% Also try different order of search coefficients.

cat_name = 'C20110808_R03';

% Sampling rate [Hz]
sampling_rate = 10000; 
 
% Load top 10 channels information.
load('Harmonics_IMs_top10.mat');

% Bipolar channel (Note this ID is not actual bipolar channel ID.)
bp_ch_id = 10;

% Frequencies of interest (harmonics and intermodulaitons)
foi_f1_harm = 46:23:250;
foi_inter = sort([177:-23:0 200+23:23:250]);
foi_now = sort([foi_f1_harm, foi_inter]);

%% Rectification 
% Load Rect data
Rect_data = load('Results_Rect_top10.mat');

% Get error from S1
S1_channel_data = Rect_data.each_area(1).each_channel;
freqs = Rect_data.each_area(1).freqs;

% Find sessoin of our interest.
search_array = repmat(cat_name, size(S1_channel_data, 2), 1);
S1_session_array = data_top10_struct(1).session;
session_this_cat = prod(S1_session_array == search_array, 2);
session_ids_this_cat = find(session_this_cat == 1);
S1_channel_data_this_cat = S1_channel_data(1, session_ids_this_cat);
S1_bipolar_chs_this_cat = data_top10_struct(1).bipolar_ch(session_ids_this_cat);

disp('Rect');
disp(cat_name);
disp(['bipolar channel: ', num2str(S1_bipolar_chs_this_cat(bp_ch_id))]);

% Get logSNR data searched with the previous method.
% Initial coeffieicent = 0. Search range from 0 to 1. Search coefficient
% alphabetically.
observed_logsnrs_now = S1_channel_data_this_cat(bp_ch_id).logsnr;
model_1_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(1).model_logsnr;
model_2_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(2).model_logsnr;
model_3_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(3).model_logsnr;
model_4_logSNR = S1_channel_data_this_cat(bp_ch_id).each_model(4).model_logsnr;

% Try other methods.
% 1. Initial coefficient = -1, -0.5, 0, 0.5, 1.0. 
% Search range from -1 to 1. Search coefficient alphabetically.

% Parameter for line search
tol = 1e-6;  
iterations = [40, 60, 60, 80]; % iterate 10 times for each coefficient
varmin = -1;
varmax = 1;

for model_id = 1:4   
    for initial_value = [-1, -0.5, 0, 0.5, 1]
        switch model_id 
            case 1 % aRect(X) + bRect(Y)
                iteration = iterations(model_id);
                model = @(a, b)model_plus(a, b);
                order = [1, 2];
                % Line search
                [L_min_var1, L_min_var2, L_error, L_errors, L_coordinate] = ...
                    fit_model_param2(model, varmin, varmax, observed_logsnrs_now, ...
                    sampling_rate, iteration, tol, initial_value, order);
                best_parameters = [L_min_var1, L_min_var2];
                % Get best signal
                model_signal = model(L_min_var1, L_min_var2);

            case 2 % aRect(X) + bRect(Y) + cRect(X)Rect(Y)
                iteration = iterations(model_id);
                model = @(a, b, c)model_without_RXY(a, b, c); 
                order = [1, 2, 3];
                % Line search
                [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                    fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                    sampling_rate, iteration, tol, initial_value, order);
                best_parameters = [L_min_var1, L_min_var2, L_min_var3];
                % Get best signal
                model_signal = model(L_min_var1, L_min_var2, L_min_var3);

            case 3 % aRect(X) + bRect(Y) + cRect(XY)
                iteration = iterations(model_id);
                model = @(a, b, c)model_without_RXRY(a, b, c);
                order = [1, 2, 3];
                % Line search
                [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                    fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                    sampling_rate, iteration, tol, initial_value, order);
                best_parameters = [L_min_var1, L_min_var2, L_min_var3];
                % Get best signal
                model_signal = model(L_min_var1, L_min_var2, L_min_var3);

            case 4 % aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY) 
                iteration = iterations(model_id);
                model = @(a, b, c, d)model_all(a, b, c, d); 
                order = [1, 2, 3, 4];
                % Line search
                [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                    fit_model_param4(model, varmin, varmax, observed_logsnrs_now, ...
                    sampling_rate, iteration, tol, initial_value, order);
                best_parameters = [L_min_var1, L_min_var2, ...
                                  L_min_var3, L_min_var4];
                % Get best signal
                model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
        end % switch model_id 

        % Compute logSNR for the best model
        [model_logsnrs, powers, freqs_] = compute_logsnrs_y(model_signal, sampling_rate);
        
        if exist('L_error_matrix') == 0
            L_error_matrix = L_error;
            best_parameters_matrix = best_parameters;
            model_signal_matrix = model_signal;
            model_logsnrs_matrix = model_logsnrs;        
        else
            L_error_matrix =[L_error_matrix; L_error];
            best_parameters_matrix = [best_parameters_matrix; best_parameters];
            model_signal_matrix = [model_signal_matrix; model_signal];
            model_logsnrs_matrix = [model_logsnrs_matrix; model_logsnrs];
        end % exist ('L_error_matrix')   
        
    end % initial_value = [-1, -0.5, 0, 0.5, 1]
    each_model_best(model_id).error = L_error_matrix;
    each_model_best(model_id).parameters = best_parameters_matrix;
    each_model_best(model_id).model_signal = model_signal_matrix;
    each_model_best(model_id).model_logsnr = model_logsnrs_matrix;
    
    clear('L_error_matrix');
    
end % model_id = 1:4

each_model_best_vary_initial = each_model_best;

% 2. Initial coefficient = 0
% Search range from 0 to 1. Search not coefficient alphabetically.

tol = 1e-6;  
varmin = -1;
varmax = 1;
iterations = [40, 60, 60, 80]; % iterate 10 times for each coefficient
initial_value = 0 ;

for model_id = 1:4        
    switch model_id 
        case 1 % aRect(X) + bRect(Y)
            iteration = iterations(model_id);
            num_order = 2;
            L_error_matrix = ones(num_order, 1);
            best_param_matrix = ones(num_order, 2);
            model = @(a, b)model_plus(a, b);
            
            % Order 1
            order_id = 1;
            order = [1, 2];
            % Line search
            [L_min_var1, L_min_var2, L_error, L_errors, L_coordinate] = ...
                fit_model_param2(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2);
            model_signal_matrix = model_signal;

            % Order 2
            order_id = 2;
            order = [2, 1];
            % Line search
            [L_min_var1, L_min_var2, L_error, L_errors, L_coordinate] = ...
                fit_model_param2(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2);
            model_signal_matrix = [model_signal_matrix; model_signal];

        case 2 % aRect(X) + bRect(Y) + cRect(X)Rect(Y)
            iteration = iterations(model_id);
            num_order = 6;
            L_error_matrix = ones(num_order, 1);
            best_param_matrix = ones(num_order, 3);
            model = @(a, b, c)model_without_RXY(a, b, c); 
            
            % Order 1
            order_id = 1;
            order = [1, 2, 3];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = model_signal;

            % Order 2            
            order_id = 2;
            order = [1, 3, 2];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 3            
            order_id = 3;
            order = [2, 1, 3];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 4            
            order_id = 4;
            order = [2, 3, 1];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;     
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];
            
            % Order 5            
            order_id = 5;
            order = [3, 1, 2];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 6            
            order_id = 6;
            order = [3, 2, 1];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];            
            
        case 3 % aRect(X) + bRect(Y) + cRect(XY)
            iteration = iterations(model_id);            
            num_order = 6;
            L_error_matrix = ones(num_order, 1);            
            best_param_matrix = ones(num_order, 3);
            model = @(a, b, c)model_without_RXRY(a, b, c);
            
            % Order 1
            order_id = 1;
            order = [1, 2, 3];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;            
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = model_signal;

            % Order 2            
            order_id = 2;
            order = [1, 3, 2];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 3            
            order_id = 3;
            order = [2, 1, 3];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 4            
            order_id = 4;
            order = [2, 3, 1];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];
            
            % Order 5            
            order_id = 5;
            order = [3, 1, 2];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 6            
            order_id = 6;
            order = [3, 2, 1];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_error, L_errors, L_coordinate] = ...
                fit_model_param3(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id,1 ) = L_error;
            best_parameters = [L_min_var1, L_min_var2, L_min_var3];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3);
            model_signal_matrix = [model_signal_matrix; model_signal];
            
        case 4 % aRect(X) + bRect(Y) + cRect(X)Rect(Y) + dRect(XY) 
            iteration = iterations(model_id);
            num_order = 4;            
            L_error_matrix = ones(num_order, 1);            
            best_param_matrix = ones(num_order, 4);
            model = @(a, b, c, d)model_all(a, b, c, d); 
            
            % Try some of them, not all.
            % Order 1
            order_id = 1;
            order = [1, 2, 3, 4];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                fit_model_param4(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2, ...
                              L_min_var3, L_min_var4];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
            model_signal_matrix = model_signal;

            % Order 2
            order_id = 2;
            order = [2, 1, 4, 3];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                fit_model_param4(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
                        
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2, ...
                              L_min_var3, L_min_var4];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
            model_signal_matrix = [model_signal_matrix; model_signal];
            
            % Order 3
            order_id = 3;
            order = [3, 2, 1, 4];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                fit_model_param4(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2, ...
                              L_min_var3, L_min_var4];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
            model_signal_matrix = [model_signal_matrix; model_signal];

            % Order 4
            order_id = 4;
            order = [4, 3, 2, 1];
            % Line search
            [L_min_var1, L_min_var2, L_min_var3, L_min_var4, L_error, L_errors, L_coordinate] = ...
                fit_model_param4(model, varmin, varmax, observed_logsnrs_now, ...
                sampling_rate, iteration, tol, initial_value, order);
            
            L_error_matrix(order_id, 1) = L_error;
            best_parameters = [L_min_var1, L_min_var2, ...
                              L_min_var3, L_min_var4];
            best_param_matrix(order_id, :) = best_parameters; 
            % Get best signal
            model_signal = model(L_min_var1, L_min_var2, L_min_var3, L_min_var4);
            model_signal_matrix = [model_signal_matrix; model_signal];
            
    end % switch model_id 
    
    for signal_id = 1:size(model_signal_matrix, 1)
        model_signal_now = model_signal_matrix(signal_id, :);
        % Compute logSNR for the best model
        [model_logsnrs, powers, freqs_] = compute_logsnrs_y(model_signal_now, sampling_rate);
        if signal_id == 1
            model_logsnrs_matrix = model_logsnrs; 
        else
            model_logsnrs_matrix = [model_logsnrs_matrix; model_logsnrs];
        end
    end % signal_id

    each_model_best(model_id).error = L_error_matrix;
    each_model_best(model_id).parameters = best_param_matrix;
    each_model_best(model_id).model_signal = model_signal_matrix;
    each_model_best(model_id).model_logsnr = model_logsnrs_matrix;

end % model_id = 1:4

each_model_best_shuffle = each_model_best;


%% Plot different versions
for model_id = 1:4
    disp(model_id);
    % Different initial coefficient
    each_model_now = each_model_best_vary_initial(model_id);
    num_versions = size(each_model_now.error, 1);
    table_now = table(each_model_now.error, each_model_now.parameters);
    table_now.Properties.VariableNames = {'best error', 'best parameters'};
    disp(table_now);
    %disp('parameters');
    %disp(each_model_now.parameters);
    %disp('best error');
    %disp(each_model_now.error);
    
    figure();
    hold on;    
    plot(freqs_, observed_logsnrs_now, 'k', 'linewidth', 3);
    for plot_id = 1:num_versions
        plot(freqs_, each_model_now.model_logsnr(plot_id, :));
    end % plot_id = 1:num_versions
    
    % Shaffle search order
    each_model_now_s = each_model_best_shuffle(model_id);
    num_versions = size(each_model_now_s.error, 1);
    table_now = table(each_model_now_s.error, each_model_now_s.parameters);
    table_now.Properties.VariableNames = {'best error', 'best parameters'};
    disp(table_now);
    %disp('parameters');
    %disp(each_model_now.parameters);
    %disp('best error');
    %disp(each_model_now.error);
            
    for plot_id = 1:num_versions
        plot(freqs_, each_model_now_s.model_logsnr(plot_id, :));
    end % plot_id = 1:num_versions
    legend();
    hold off   
    xlim([0, 250]);
    xlabel('freqency [Hz]');
    ylabel('logSNR [dB]');
    
end % model_id = 1:4

%% Functions
function [min_var1, min_var2, error, errors, coordinate] = fit_model_param2(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol, initial_value, order)
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
%           initial_value = initial coefficient before searching it. 
%           order = list of numbers, indicating the order of searching 
%                   coefficient. a=1 and b=2
%                   e.g. list = [1, 2] -> search a first then b.
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
            target = order(1); %  
            if i_iter == 1
                varfixed = initial_value; % 
            else 
                varfixed = coordinate(i_iter - 1, order(2));
            end % i_iter == 1
            
            [new_var, error]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        elseif mod(i_iter, 2) == 0 % Then fit var2 
            target = order(2);
            varfixed = coordinate(i_iter - 1, order(1));
            [new_var, error]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
        end % mod(i_iter, 2) == 1  
        
        switch target
            case 1
                coordinate_now = [new_var, varfixed]; % new coordinate 
            case 2 
                coordinate_now = [varfixed, new_var];
        end

        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = error;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
end

function [min_var1, min_var2, min_var3, error, errors, coordinate] = fit_model_param3(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol, initial_value, order)
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
%           initial_value = initial coefficient before searching it. 
%           order = list of numbers, indicating the order of searching 
%                   coefficient. a=1, b=2, and c=3
%                   e.g. list = [1, 2, 3] -> search a first, then b,
%                   finally c. 
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
            target = order(1);
            if i_iter == 1
                varfixed = [initial_value, initial_value]; 
            else 
                switch target
                    case 1
                        varfixed = coordinate(i_iter - 1, 2:3);
                    case 2
                        varfixed = [coordinate(i_iter - 1, 1),...
                                    coordinate(i_iter - 1, 3)];
                    case 3
                        varfixed = coordinate(i_iter - 1, 1:2);                        
                end % switch target                              
            end % i_iter == 1
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        elseif mod(i_iter, 3) == 2 % Then fit var2 
            target = order(2);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:3);
                case 2
                    varfixed = [coordinate(i_iter - 1, 1),...
                                coordinate(i_iter - 1, 3)];
                case 3
                    varfixed = coordinate(i_iter - 1, 1:2);                        
            end % switch target                              
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
        elseif mod(i_iter, 3) == 0
            target = order(3);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:3);
                case 2
                    varfixed = [coordinate(i_iter - 1, 1),...
                                coordinate(i_iter - 1, 3)];
                case 3
                    varfixed = coordinate(i_iter - 1, 1:2);                        
            end % switch target                              
            [new_var, error]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        end % mod(i_iter, 2) == 1  
        
        % Update corrdinate
        switch target
            case 1
                coordinate_now = [new_var, varfixed(1), varfixed(2)]; % new coordinate 
            case 2 
                coordinate_now = [varfixed(1), new_var, varfixed(2)]; % new coordinate                 
            case 3 
                coordinate_now = [varfixed(1), varfixed(2), new_var]; % new coordinate                 
        end % switch target

        
        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = error;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
    min_var3 = coordinate_now(3);
end

function [min_var1, min_var2, min_var3, min_var4, error, errors, coordinate] = fit_model_param4(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, tol, initial_value, order)
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
%           initial_value = initial coefficient before searching it. 
%           order = list of numbers, indicating the order of searching 
%                   coefficient. a=1, b=2, c=3, and d=4
%                   e.g. list = [1, 2, 3, 4] -> search a first, then b,
%                   then c, finally d. 
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
            target = order(1);
            if i_iter == 1
                varfixed = [initial_value, initial_value, initial_value]; 
            else
                switch target
                    case 1
                        varfixed = coordinate(i_iter - 1, 2:4);
                    case 2
                        varfixed = [coordinate(i_iter - 1, 1), ...
                                    coordinate(i_iter - 1, 3:4)];
                    case 3
                        varfixed = [coordinate(i_iter - 1, 1:2), ...
                                    coordinate(i_iter - 1, 4)];                    
                    case 4
                        varfixed = coordinate(i_iter - 1, 1:3);
                end % switch target                            
            end % i_iter == 1
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        elseif mod(i_iter, 4) == 2 % Fit var2 
            target = order(2);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:4);
                case 2
                    varfixed = [coordinate(i_iter - 1, 1), ...
                                coordinate(i_iter - 1, 3:4)];
                case 3
                    varfixed = [coordinate(i_iter - 1, 1:2), ...
                                coordinate(i_iter - 1, 4)];                    
                case 4
                    varfixed = coordinate(i_iter - 1, 1:3);
            end % switch target            
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);       
        elseif mod(i_iter, 4) == 3 % Fit var3
            target = order(3);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:4);
                case 2
                    varfixed = [coordinate(i_iter - 1, 1), ...
                                coordinate(i_iter - 1, 3:4)];
                case 3
                    varfixed = [coordinate(i_iter - 1, 1:2), ...
                                coordinate(i_iter - 1, 4)];                    
                case 4
                    varfixed = coordinate(i_iter - 1, 1:3);
            end % switch target
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        elseif mod(i_iter, 4) == 0
            target = order(4);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:4);
                case 2
                    varfixed = [coordinate(i_iter - 1, 1), ...
                                coordinate(i_iter - 1, 3:4)];
                case 3
                    varfixed = [coordinate(i_iter - 1, 1:2), ...
                                coordinate(i_iter - 1, 4)];                    
                case 4
                    varfixed = coordinate(i_iter - 1, 1:3);
            end % switch target
            [new_var, error]= golden_section_parm4(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol);
        end % mod(i_iter, 4) == 1  
        
        % Update coordinate
        switch target 
            case 1
                coordinate_now = [new_var, varfixed(1), varfixed(2), varfixed(3)]; % new coordinate 
            case 2
                coordinate_now = [varfixed(1), new_var, varfixed(2), varfixed(3)]; % new coordinate            
            case 3
                coordinate_now = [varfixed(1), varfixed(2), new_var, varfixed(3)]; % new coordinate                         
            case 4
                coordinate_now = [varfixed(1), varfixed(2), varfixed(3), new_var]; % new coordinate                                     
        end
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

