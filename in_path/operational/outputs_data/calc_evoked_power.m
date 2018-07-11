function data_out = calc_evoked_power(data, baseline, flag)

% calc_evoked_power: Calculates the evoked power
%   The evoked power is defined as the amount of power for a given channel
%   less the mean of the null condition. 
%   
%   data_out = calc_evoked_power(data, baseline, flag) calculates the
%   evoked power contained in data.trial. flag determines whether it
%   calcuates the baseline or other condition. Set flag = 1 if it is the
%   null/baseline condition, or set to 0 to run other conditions. Default
%   is 0 (when there are only 2 arguments).
%   baseline is ignored when flag == 1. Recommendation: set to [].

% Written by Rannee Li, Jun 2017.


% set default flag behaviour
if nargin<2
    flag = 0;
end

% select from the subfunctions

if flag
    data_out = cep_base(data);
else
    data_out = cep_nonbase(data, baseline);
end

% update filename
data_out.custom.filename = [data.custom.filename '_evkdpwr'];

end
% ----- SUBFUNCTIONS BELOW -----

function data_out = cep_base(data)
data_out = data;

% calculate the mean of the baseline across trials
trialmean = mean(data.trial, 3);

% repmat baseline
data_out.trial = repmat(trialmean, [1, 1, data.custom.ntrials]);

end


function data_out = cep_nonbase(data, baseline)
data_out = data;

% subtract baseline
data_out.trial = data.trial - baseline.trial;

end