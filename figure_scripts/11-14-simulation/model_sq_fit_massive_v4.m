function model_sq_fit_massive_v4(channel_id, area_id)
%%Fit all Sq models to a channel data. 
% Error = Sum across all frequencies.
% Input:    channel_id = bipolar channel id. (each area)
%           area_id = 1 or 2 (i.e. S1 or S2)
% Output:   Save fitting results.


%% Parameters
% channel_id = 99;
% area_id = 1;

% Sampling rate [Hz]
sampling_rate = 10000; 

% Add path
addpath(genpath('../../in_path'));
addpath(genpath('../../scripts'));

%% Load top 10 channels information.
loaded_data = load('Harmonics_IMs_top10.mat');

% Load data for this channel in this area 
observed_logsnr = loaded_data.data_top10_struct(area_id).mean_logsnrs(channel_id, :);
session = loaded_data.data_top10_struct(area_id).session(channel_id, :);
bipolar_ch = loaded_data.data_top10_struct(area_id).bipolar_ch(channel_id);

%% Rectification fittings
% Get logSNR data searched with the previous method.
% Initial coefficient = [-1, -0.5, 0, 0.5, 1];
% Initial range range from -1 to 1. 
% Search all order 
varmin = -1;
varmax = 1;
num_coefficients = [2, 3];
iterations = num_coefficients*15 ; % iterate 5 times for each coefficient
initial_value_candidates = [-1, 0, 1];

tic; % how log to compute all models for one channel
for model_id = 1:2
    id_incre = 0;   
    iteration = iterations(model_id); 
    initial_values_matrix = permn(initial_value_candidates, ...
                                  num_coefficients(model_id));
    for init_ind = 1:size(initial_values_matrix, 1)
        initial_values = initial_values_matrix(init_ind, :);
        switch model_id 
            case 1 % aSq(X) + bSq(Y)
                order_matrix = perms([2, 1]);
                num_order = size(order_matrix, 1);
                model = @(a, b)model_plus_sq(a, b);
                
                for order_id = 1:num_order
                    order = order_matrix(order_id, :);                    

                    % Line search
                    [L_min_var1, L_min_var2, L_error, ~, ~] = ...
                        fit_model_param2(model, varmin, varmax, observed_logsnr, ...
                        sampling_rate, iteration, initial_values, order);
                    best_parameters = [L_min_var1, L_min_var2];
                    
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2);
                    
                    % Keep data
                    if init_ind == 1 && order_id == 1
                        L_error_matrix = nan(num_order*size(initial_values_matrix, 1), 1);
                        best_param_matrix = nan(num_order*size(initial_values_matrix, 1), 2);
                        model_signal_matrix = nan(num_order*size(initial_values_matrix, 1), ...
                                                    length(model_signal));
                    end % init_ind == 1 && order_id == 1
                    id_incre = id_incre + 1;
                    L_error_matrix(id_incre, 1) = L_error;
                    best_param_matrix(id_incre, :) = best_parameters; 
                    model_signal_matrix(id_incre, :) = model_signal;
                end % order_id = 1:num_order
                              
            case 2 % aSq(X) + bSq(Y) + cSq(X)Sq(Y)
                order_matrix = perms([3, 2, 1]);
                num_order = size(order_matrix, 1);
                model = @(a, b, c)model_without_RXY_sq(a, b, c); 

                for order_id = 1:num_order
                    order = order_matrix(order_id, :);
                    
                    % Line search
                    [L_min_var1, L_min_var2, L_min_var3, L_error, ~, ~] = ...
                        fit_model_param3(model, varmin, varmax, observed_logsnr, ...
                        sampling_rate, iteration, initial_values, order);
                    best_parameters = [L_min_var1, L_min_var2, L_min_var3];
                    
                    % Get best signal
                    model_signal = model(L_min_var1, L_min_var2, L_min_var3);
                    
                    % Keep data
                    if init_ind == 1 && order_id == 1
                        L_error_matrix = nan(num_order*size(initial_values_matrix, 1), 1);
                        best_param_matrix = nan(num_order*size(initial_values_matrix, 1), 3);                        
                        model_signal_matrix = nan(num_order*size(initial_values_matrix, 1), ...
                                                    length(model_signal));
                    end % init_ind == 1 && order_id == 1      
                    id_incre = id_incre + 1;
                    L_error_matrix(id_incre,1 ) = L_error;            
                    best_param_matrix(id_incre, :) = best_parameters;                     
                    model_signal_matrix(id_incre, :) = model_signal;
                end % order_id = 1:num_order                                               
        end % switch model_id 

        for signal_id = 1:size(model_signal_matrix, 1)
            model_signal_now = model_signal_matrix(signal_id, :);
            % Compute logSNR for the best model
            [model_logsnrs, ~, freqs] = compute_logsnrs_y(model_signal_now, sampling_rate);
            if signal_id == 1
                model_logsnrs_matrix = model_logsnrs; 
            else
                model_logsnrs_matrix = [model_logsnrs_matrix; model_logsnrs];
            end
        end % signal_id
        
    end % initial_value
    % Store results from all settings 
    each_model_all(model_id).error = L_error_matrix;
    each_model_all(model_id).parameters = best_param_matrix;
    each_model_all(model_id).model_signal = model_signal_matrix;
    each_model_all(model_id).model_logsnr = model_logsnrs_matrix;

    % Store the best across all settings
    [~, min_ind] = min(L_error_matrix);
    each_model_best(model_id).error = L_error_matrix(min_ind);
    each_model_best(model_id).parameters = best_param_matrix(min_ind,:);
    each_model_best(model_id).model_signal = model_signal_matrix(min_ind, :);
    each_model_best(model_id).model_logsnr = model_logsnrs_matrix(min_ind, :);   
               
end % model_id = 1:4
toc; % % how log to compute all models for one channel

% Save data 
% all settings 
dir_name_all_settings = './sq_top10_results_all_settings_v4/';
if ~exist(dir_name_all_settings, 'dir')
    mkdir(dir_name_all_settings);
end % ~exist(dir_name_all_settings, 'dir');

file_name_all_settings = ['all_settings_sq_S', num2str(area_id), '_id_', ...
             num2str(channel_id), '.mat'];
save([dir_name_all_settings, file_name_all_settings], 'each_model_all', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');

% only the best setting
dir_name_best_setting = './sq_top10_results_best_setting_v4/';
if ~exist(dir_name_best_setting, 'dir')
    mkdir(dir_name_best_setting);
end % ~exist(dir_name_all_settings, 'dir');
file_name_best_setting = ['best_setting_sq_S', num2str(area_id), '_id_', ...
             num2str(channel_id), '.mat'];
save([dir_name_best_setting, file_name_best_setting], 'each_model_best', ...
    'observed_logsnr', 'freqs', 'area_id', 'session', 'bipolar_ch');

%% Functions
function [min_var1, min_var2, err, errors, coordinate] = fit_model_param2(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, initial_values, order)
% Find var1 and var 2 giving the local minimum by line serach. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  func = function to be accessed. Consider two-variable
%                  function (Model)
%           varmin = min of variable 
%           varmax = max of variable 
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
%           iteration = how many times you want to repeat linear search
%           initial_values = initial coefficient before searching them. 
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
    tols = [1, 1e-1, 1e-2];
    increment = 0;
    
    for i_iter = 1:iteration
        if increment < 5
            tol = tols(1);
        elseif increment < 10
            tol = tols(2);
        else 
            tol = tols(3);
        end % increment 
            
        if mod(i_iter, 2) == 1 % Fit var1 first
            increment = increment + 1;
            target = order(1); %  
            if i_iter == 1
                switch target 
                    case 1
                        varfixed = initial_values(2); % 
                        varmin_now = varmin + initial_values(1);
                        varmax_now = varmax + initial_values(1);
                    case 2
                        varfixed = initial_values(1); % 
                        varmin_now = varmin + initial_values(2);
                        varmax_now = varmax + initial_values(2);
                end % switch target          
            else
                varfixed = coordinate(i_iter - 1, order(2));
                varmin_now = coordinate(i_iter - 1, order(1));
                varmax_now = varmax + coordinate(i_iter - 1, order(1));                
            end % i_iter == 1
            [new_var, err]= golden_section_parm2(model, datalogsnrs, ...
                sampling_rate, varmin_now, varmax_now, varfixed, target, tol);
        elseif mod(i_iter, 2) == 0 % Then fit var2 
            target = order(2);
            varfixed = coordinate(i_iter - 1, order(1));
            varmin_now = coordinate(i_iter - 1, order(2));
            varmax_now = varmax + coordinate(i_iter - 1, order(2));
            
            [new_var, err]= golden_section_parm2(model, datalogsnrs, ...
                sampling_rate, varmin_now, varmax_now, varfixed, target, tol);       
        end % mod(i_iter, 2) == 1  
        
        switch target
            case 1
                coordinate_now = [new_var, varfixed]; % new coordinate 
            case 2 
                coordinate_now = [varfixed, new_var];
        end

        % Store data
        coordinate(i_iter, :) = coordinate_now;
        errors(i_iter, 1) = err;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
end

function [min_var1, min_var2, min_var3, err, errors, coordinate] = fit_model_param3(model, varmin, varmax, datalogsnrs, sampling_rate, iteration, initial_values, order)
% Find var1 and var 2 giving the local minimum by line serach. We
% minimise error (difference of logSNR between model and data across foi).
%   Input:  func = function to be accessed. Consider two-variable
%                  function (Model)
%           varmin = min of variable 
%           varmax = max of variable 
%           datalogsnrs = data that we fit the function to.
%           sampling_rate = sampling rate for func
%           iteration = how many times you want to repeat linear search
%           initial_values = initial coefficient before searching them. 
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
    tols = [1, 1e-1, 1e-2];
    increment = 0;
    
    for i_iter = 1:iteration
        if increment < 5
            tol = tols(1);
        elseif increment < 10 
            tol = tols(2);
        else 
            tol = tols(3);
        end % increment 

        if mod(i_iter, 3) == 1 % Fit var1 first
            increment = increment + 1;
            target = order(1);
            if i_iter == 1
                switch target 
                    case 1
                        varfixed = initial_values(2:3); 
                        varmin_now = varmin + initial_values(1);
                        varmax_now = varmax + initial_values(1);
                    case 2
                        varfixed = [initial_values(1),...
                                    initial_values(3)]; 
                        varmin_now = varmin + initial_values(2);
                        varmax_now = varmax + initial_values(2);
                    case 3
                        varfixed = initial_values(1:2);  
                        varmin_now = varmin + initial_values(3);
                        varmax_now = varmax + initial_values(3);
                end % switch target                       
            else 
                switch target
                    case 1
                        varfixed = coordinate(i_iter - 1, 2:3);
                        varmin_now = coordinate(i_iter - 1, 1);
                        varmax_now = varmax + coordinate(i_iter - 1, 1);                
                    case 2
                        varfixed = [coordinate(i_iter - 1, 1),...
                                    coordinate(i_iter - 1, 3)];                                
                        varmin_now = coordinate(i_iter - 1, 2);
                        varmax_now = varmax + coordinate(i_iter - 1, 2);
                    case 3
                        varfixed = coordinate(i_iter - 1, 1:2);                         
                        varmin_now = coordinate(i_iter - 1, 3);
                        varmax_now = varmax + coordinate(i_iter - 1, 3);
                end % switch target                              
            end % i_iter == 1
            [new_var, err]= golden_section_parm3(model, datalogsnrs, ...
                sampling_rate, varmin_now, varmax_now, varfixed, target, tol);
        elseif mod(i_iter, 3) == 2 % Then fit var2 
            target = order(2);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:3);
                    varmin_now = coordinate(i_iter - 1, 1);
                    varmax_now = varmax + coordinate(i_iter - 1, 1);                
                case 2
                    varfixed = [coordinate(i_iter - 1, 1),...
                                coordinate(i_iter - 1, 3)];                                
                    varmin_now = coordinate(i_iter - 1, 2);
                    varmax_now = varmax + coordinate(i_iter - 1, 2);
                case 3
                    varfixed = coordinate(i_iter - 1, 1:2);                         
                    varmin_now = coordinate(i_iter - 1, 3);
                    varmax_now = varmax + coordinate(i_iter - 1, 3);
            end % switch target                              
            [new_var, err]= golden_section_parm3(model, datalogsnrs, ...
                sampling_rate, varmin_now, varmax_now, varfixed, target, tol);       
        elseif mod(i_iter, 3) == 0
            target = order(3);
            switch target
                case 1
                    varfixed = coordinate(i_iter - 1, 2:3);
                    varmin_now = coordinate(i_iter - 1, 1);
                    varmax_now = varmax + coordinate(i_iter - 1, 1);                
                case 2
                    varfixed = [coordinate(i_iter - 1, 1),...
                                coordinate(i_iter - 1, 3)];                                
                    varmin_now = coordinate(i_iter - 1, 2);
                    varmax_now = varmax + coordinate(i_iter - 1, 2);
                case 3
                    varfixed = coordinate(i_iter - 1, 1:2);                         
                    varmin_now = coordinate(i_iter - 1, 3);
                    varmax_now = varmax + coordinate(i_iter - 1, 3);
            end % switch target                              
            [new_var, err]= golden_section_parm3(model, datalogsnrs, ...
                sampling_rate, varmin_now, varmax_now, varfixed, target, tol);
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
        errors(i_iter, 1) = err;
    end
    % Return the last one
    min_var1 = coordinate_now(1);
    min_var2 = coordinate_now(2);    
    min_var3 = coordinate_now(3);
end

function [var, err, xs, ys]= golden_section_parm2(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
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
%           xs = history of x (search for a coefficient)
%           ys = history of y (error)
%   Note:  Data should have the same sampling rate with func.
%   
    r = (sqrt(5) - 1)/2;

    % Set initial a,b,c,d and their func(x) values i.e. error(x)
    a = varmin;
    d = varmax;
    b = r*a + (1-r)*d;  % 1-r:r
    c = (1-r)*a + r*d;  % r:1-r
    
    if target == 1 % find var1 and fix var2
        signal_a = model(a, varfixed(1));
        signal_b = model(b, varfixed(1));
        signal_c = model(c, varfixed(1));
        signal_d = model(d, varfixed(1));
    elseif target == 2 % fix var1 and find var2
        signal_a = model(varfixed(1), a);
        signal_b = model(varfixed(1), b);
        signal_c = model(varfixed(1), c);
        signal_d = model(varfixed(1), d);
    end
    f_a = compute_error(signal_a, datalogsnrs, sampling_rate);
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
    f_d = compute_error(signal_d, datalogsnrs, sampling_rate);
    
    % threshold to stop while loop
    % Do at least one loop 
    tol_ = inf; 
    
    xs = [a, b, c, d];
    y_list = [f_a, f_b, f_c, f_d];    
    ys = y_list;
    
    % Update variabel until (d - a) < tol
    while (0.001 < (d-a))||(tol < tol_)
        [~, min_ind] = min(y_list);
        switch min_ind
            case 1 % Extend search range toward smaller area.
                % Keep the previous data temporarily
                a_ = a; b_ = b; d_ = d; %c_ = c;
                f_a_ = f_a; f_b_ = f_b; f_d_ = f_d; %f_c_ = f_c;
                % Update coefficients
                a = a_ - (d_ - b_);
                if target == 1 % update var1 and fix var2 
                    signal_a = model(a,  varfixed(1));
                elseif target == 2 % fix var1 and update var2
                    signal_a = model(varfixed(1), a);
                end
                f_a = compute_error(signal_a, datalogsnrs, sampling_rate);
                b = a_;
                f_b = f_a_;
                c = b_; 
                f_c = f_b_; 
                d = d_;
                f_d = f_d_;
            case 2 % Eliminate x > c by replacing d with c.
                % Keep the previous data temporarily
                a_ = a; b_ = b; c_ = c; %d_ = d;
                f_a_ = f_a; f_b_ = f_b; f_c_ = f_c; %f_d_ = f_d;
                % Update coefficients
                a = a_;
                f_a = f_a_;
                b = r*a_ + (1-r)*c_;  % 1-r:r
                if target == 1 % update var1 and fix var2 
                    signal_b = model(b, varfixed(1));  
                elseif target == 2 % fix var1 and update var2
                    signal_b = model(varfixed(1), b);  
                end
                f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
                c = b_; 
                f_c = f_b_;
                d = c_;
                f_d = f_c_;
            case 3 % Eliminate x < b by replacing a with b.
                % Keep the previous data temporarily
                b_ = b; c_ = c; d_ = d; %a_ = a; 
                f_b_ = f_b; f_c_ = f_c; f_d_ = f_d; %f_a_ = f_a; 
                % Update coefficients
                a = b_;
                f_a = f_b_;
                b = c_;
                f_b = f_c_;
                c = (1-r)*b_ + r*d_;  % r:1-r
                if target == 1 % update var1 and fix var2 
                    signal_c = model(c, varfixed(1));  
                elseif target == 2 % fix var1 and update var2
                    signal_c = model(varfixed(1), c);  
                end           
                f_c = compute_error(signal_c, datalogsnrs, sampling_rate);                
                d = d_;
                f_d = f_d_;
            case 4 % Extend search range toward smaller area.
                % Keep the previous data temporarily
                a_ = a; c_ = c; d_ = d; %b_ = b; 
                f_a_ = f_a; f_c_ = f_c; f_d_ = f_d; %f_b_ = f_b; 
                % Update coefficients
                a = a_;
                f_a = f_a_;
                b = c_;
                f_b = f_c_;
                c = d_;
                f_c = f_d_;
                d = d_ + (c_ - a_);
                if target == 1 % update var1 and fix var2 
                    signal_d = model(d, varfixed(1));
                elseif target == 2 % fix var1 and update var2
                    signal_d = model(varfixed(1), d);
                end
                f_d = compute_error(signal_d, datalogsnrs, sampling_rate);
        end % switch min_ind
        
        % Keep search history
        x_list = [a, b, c, d];
        y_list = [f_a, f_b, f_c, f_d];
                
        % threshold to stop while loop
        tol_ = min(ys(end,:)) - min(y_list);
        
        xs = [xs; x_list];
        ys = [ys; y_list];
        if a == b || c == d
            break;
        end
    end % while 
    
    % Get the best var and its error
    [~, min_ind] = min(y_list);
    var = x_list(min_ind);
    err = y_list(min_ind);
 end

function [var, err, xs, ys]= golden_section_parm3(model, datalogsnrs, sampling_rate, varmin, varmax, varfixed, target, tol)
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
%           xs = history of x (search for a coefficient)
%           ys = history of y (error)

%   Note:  Data should have the same sampling rate with func.
%   
    r = (sqrt(5) - 1)/2;

    % Set initial a,b,c,d and their func(x) values i.e. error(x)
    a = varmin;
    d = varmax;
    b = r*a + (1-r)*d;  % 1-r:r
    c = (1-r)*a + r*d;  % r:1-r
    
    if target == 1 % find var1 and fix var2 and var3
        signal_a = model(a, varfixed(1), varfixed(2));
        signal_b = model(b, varfixed(1), varfixed(2));
        signal_c = model(c, varfixed(1), varfixed(2));
        signal_d = model(d, varfixed(1), varfixed(2));        
    elseif target == 2 % fix var1 and var3 and find var2
        signal_a = model(varfixed(1), a, varfixed(2));
        signal_b = model(varfixed(1), b, varfixed(2));
        signal_c = model(varfixed(1), c, varfixed(2));
        signal_d = model(varfixed(1), d, varfixed(2));
    elseif target == 3 % fix var1 and var2 and find var3
        signal_a = model(varfixed(1), varfixed(2), a);
        signal_b = model(varfixed(1), varfixed(2), b);
        signal_c = model(varfixed(1), varfixed(2), c);
        signal_d = model(varfixed(1), varfixed(2), d);
    end
    
    f_a = compute_error(signal_a, datalogsnrs, sampling_rate);
    f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
    f_c = compute_error(signal_c, datalogsnrs, sampling_rate);
    f_d = compute_error(signal_d, datalogsnrs, sampling_rate);
    
    % threshold to stop while loop
    % Do at least one loop 
    tol_ = inf; 
    
    xs = [a, b, c, d];
    y_list = [f_a, f_b, f_c, f_d];    
    ys = y_list;
    
    % Update variabel until (d - a) < tol
    while  (0.001 < (d-a))||(tol < tol_)
        [~, min_ind] = min(y_list);
        switch min_ind
            case 1 % Extend search range toward smaller area.
                % Keep the previous data temporarily
                a_ = a; b_ = b; d_ = d; %c_ = c; 
                f_a_ = f_a; f_b_ = f_b; f_d_ = f_d; %f_c_ = f_c; 
                % Update coefficients
                a = a_ - (d_ - b_);
                if target == 1 % update var1 and fix var2 and var3
                    signal_a = model(a, varfixed(1), varfixed(2));  
                elseif target == 2 % fix var1 and var3 and update var2
                    signal_a = model(varfixed(1), a, varfixed(2));  
                elseif target == 3 % fix var1 and var2 and update var3
                    signal_a = model(varfixed(1), varfixed(2), a);
                end
                f_a = compute_error(signal_a, datalogsnrs, sampling_rate);
                b = a_;
                f_b = f_a_;
                c = b_; 
                f_c = f_b_; 
                d = d_;
                f_d = f_d_;
            case 2 % Eliminate x > c by replacing d with c.
                % Keep the previous data temporarily
                a_ = a; b_ = b; c_ = c; %d_ = d;
                f_a_ = f_a; f_b_ = f_b; f_c_ = f_c; %f_d_ = f_d;
                % Update coefficients
                a = a_;
                f_a = f_a_;
                b = r*a_ + (1-r)*c_;  % 1-r:r
                if target == 1 % update var1 and fix var2 and var3
                    signal_b = model(b, varfixed(1), varfixed(2));  
                elseif target == 2 % fix var1 and var3 and update var2
                    signal_b = model(varfixed(1), b, varfixed(2));  
                elseif target == 3 % fix var1 and var2 and update var3
                    signal_b = model(varfixed(1), varfixed(2), b);                
                end           
                f_b = compute_error(signal_b, datalogsnrs, sampling_rate);
                c = b_; 
                f_c = f_b_;
                d = c_;
                f_d = f_c_;
            case 3 % Eliminate x < b by replacing a with b.
                % Keep the previous data temporarily
                b_ = b; c_ = c; d_ = d; %a_ = a; 
                f_b_ = f_b; f_c_ = f_c; f_d_ = f_d; %f_a_ = f_a; 
                % Update coefficients
                a = b_;
                f_a = f_b_;
                b = c_;
                f_b = f_c_;
                c = (1-r)*b_ + r*d_;  % r:1-r
                if target == 1 % update var1 and fix var2 and var3
                    signal_c = model(c, varfixed(1), varfixed(2));  
                elseif target == 2 % fix var1 and var3 and update var2
                    signal_c = model(varfixed(1), c, varfixed(2));  
                elseif target == 3 % fix var1 and var2 and update var3
                    signal_c = model(varfixed(1), varfixed(2), c);                
                end           
                f_c = compute_error(signal_c, datalogsnrs, sampling_rate);                
                d = d_;
                f_d = f_d_;
            case 4 % Extend search range toward smaller area.
                % Keep the previous data temporarily
                a_ = a; c_ = c; d_ = d; %b_ = b; 
                f_a_ = f_a; f_c_ = f_c; f_d_ = f_d; %f_b_ = f_b; 
                % Update coefficients
                a = a_;
                f_a = f_a_;
                b = c_;
                f_b = f_c_;
                c = d_;
                f_c = f_d_;
                d = d_ + (c_ - a_);
                if target == 1 % update var1 and fix var2 and var3
                    signal_d = model(d, varfixed(1), varfixed(2));  
                elseif target == 2 % fix var1 and var3 and update var2
                    signal_d = model(varfixed(1), d, varfixed(2));  
                elseif target == 3 % fix var1 and var2 and update var3
                    signal_d = model(varfixed(1), varfixed(2), d);                
                end           
                f_d = compute_error(signal_d, datalogsnrs, sampling_rate);
        end % switch min_ind
        
        % Keep search history
        x_list = [a, b, c, d];
        y_list = [f_a, f_b, f_c, f_d];
                
        % threshold to stop while loop
        tol_ = min(ys(end,:)) - min(y_list);
        
        xs = [xs; x_list];
        ys = [ys; y_list];   
        if a == b || c == d
            break;
        end
    end % while 
    
    % Get the best var and its error
    [~, min_ind] = min(y_list);
    var = x_list(min_ind);
    err = y_list(min_ind);

end

function err = compute_error(model_signal, datalogsnrs, sampling_rate)
% Compute error
%   Input: sampling_rate
%   Output: err

    % Compute logSNR
    [modellogsnrs, ~, freqs] = compute_logsnrs_y(model_signal, sampling_rate);
    
    [~, foi_ind] = find_closest(freqs, 250);
    err = sum(abs(modellogsnrs(1, 1:foi_ind) - datalogsnrs(1, 1:foi_ind)));
    
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

%% Model 
function model = model_plus_sq(a, b)
% aSq(X) + bSq(Y)
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
    
    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;
    
    % model
    model = (a*Sq_comp1) + (b*Sq_comp2);
end

function model = model_without_RXY_sq(a, b, c)
% aSq(X) + bSq(Y) + cSq(X)Sq(Y)
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

    % Square 
    Sq_comp1 = comp1;
    Sq_comp1(sin(comp1) > 0) = 1;
    Sq_comp1(sin(comp1) <=0) = -1;

    Sq_comp2 = comp2;
    Sq_comp2(sin(comp2) > 0) = 1;
    Sq_comp2(sin(comp2) <=0) = -1;
        
    % model 
    model = a*Sq_comp1 + b*Sq_comp2 + c*Sq_comp1 .* Sq_comp2;
end

end