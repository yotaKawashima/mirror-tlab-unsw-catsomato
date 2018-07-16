function data_out = resample_chooser(data, method, otherarg)
% resample_chooser: Resamples data using a selected method.
%   data_out = resample_chooser(data, 'closest', poi) resamples the data 
%   using the function 'dwnsmpl'. 
%
%   data_out = resample_chooser(data, method, target) performs the
%   preprocessing needed and calls interp1 for any of the methods that it
%   can take ('nearest', 'next', 'previous', 'linear','spline','pchip', 
%   'cubic').

% Written by Rannee Li, Feb 2015
% Update notes:
%   Jul 2015: Fixed help
%   Jun 2017: Updated the subfunction so that only one copy of time saved.

switch method
    case 'closest'
        data_out = dwnsmpl(data, otherarg);
        
    case 'resample'
        % as this is matlab-inbuilt, need to process the data before it can
        % be used.
        warning('Implementation in progress. Errors may occur.')
        
        data_out = resample(data.trial);
        
    case {'nearest', 'next', 'previous', 'linear','spline','pchip', 'cubic'}
        % use spline to interpolate. this is also inbuilt and thus requires
        % pre- and post-processing, which is done in a subfunction.
        
        data_out = rsc_inbuilthelper(data, otherarg, method);
        
        
    otherwise
        error('method argument invalid.')
end

% -------------------------- SUBFUNCTIONS FOLLOW -------------------------- 

function data_out = rsc_inbuilthelper(data, target, method)

% first, copy the data
data_out = data;

% preallocate
data_out.trial = zeros(data.custom.nsignals, numel(target), data.custom.ntrials);

times = data.time{1}; %repmat(data.time{1}, [data.custom.ntrials, 1]);

% have to loop as this isn't fully vectorised

    for t = 1:data.custom.ntrials
        for s = 1:data.custom.nsignals
        yy = interp1(times, data.trial(s, :, t), target, method);
        %yy = spline(data.time{1},data.trial(s, :, t),target);
        
        data_out.trial(s, :, t) = yy;
        end
        
    end
data_out.time{1} = target;

% fix other variables
data_out.custom.nsamples        = numel(target);
data_out.fsample                = 1/(target(2) - target(1));
data_out.custom.filename        = [data.custom.filename '_rsamp' method(1) method(3)];
data_out.custom.resamp_method   = method;