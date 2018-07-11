function [freq_keep, keep_i, hars] = remove_tags(stimf, frange, freq_vect, k, T, output_type)

% remove_tags: find indices that are needed to remove response frequencies
%




% modified from no_resp_indices.m

%% parse inputs and set defaults
% frange
if isempty(frange)
    frange = [min(freq_vect(:)), max(freq_vect(:))];
end

% freq_vect
if iscell(freq_vect)
    warning('Detected freq_vect is a cell. Using first element.')
    freq_vect = freq_vect{1};
end

% k
if nargin < 4 || isempty(k)
    k = 1;
end

% T
if nargin < 5 || isempty(T)
    T = 2;
end

% output_type
if nargin < 6
    output_type = true;
    % true for vector of indices, false for vector of logicals
elseif strcmp(output_type, 'vector')
    output_type = true;
elseif strcmp(output_type, 'logical')
    output_type = false;
else
    warning('output_type flag set to illegal value. Assuming ''vector''.')
    output_type = true;
end

%%

lind = find_closest_ind(freq_vect, frange(1));
hind = find_closest_ind(freq_vect, frange(2));

hars = [stimf(1):stimf(1):frange(2) stimf(2):-stimf(1):frange(1) stimf(2):stimf(1):frange(2)];
%       harmonics                   intermodulation              intermodulation pt 2
hars = unique(hars(hars>frange(1) & hars<frange(2)));

% remove the freq +- half the bandwidth
halfrmband = (k+1)/(2*T); 

rmf_lb = hars-halfrmband;
rmf_ub = hars+halfrmband;

keep_ind = lind:hind;

f_base = freq_vect(lind:hind);
% remove response frequencies
for m = 1:length(rmf_lb)
    keep_ind(f_base>rmf_lb(m) & f_base<rmf_ub(m)) = 0;
end

keep_i = keep_ind>0;

if output_type
    keep_i = keep_ind(keep_i);
end

freq_keep = freq_vect(keep_i);

if nargout < 3
    clear hars
end