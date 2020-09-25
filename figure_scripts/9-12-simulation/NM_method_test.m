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
    coordinate = [0, 0; v
        
    ar1max, 0; 0, var2max];
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
