function [x_mean, CI, err] = calc_ci(x, percentage, std_flag)

tail = (100-percentage)/2/100;

SEM = std(x, std_flag)/sqrt(length(x));    % Standard Error
ts = tinv([tail, 1-tail],length(x)-1);            % T-Score

x_mean = mean(x);                         % mean

err = ts*SEM;

CI = x_mean + err;                               % Confidence Intervals

err = abs(err(2));
if nargout <3
    clear err
end