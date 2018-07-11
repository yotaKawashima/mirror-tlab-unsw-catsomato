function [freq_keep, keep_i, hars] = no_resp_indices(stimf, frange, freq_vect, output_type)


if nargin < 4
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


if isempty(frange)
    frange = [min(freq_vect(:)), max(freq_vect(:))];
end


lind = find_closest_ind(freq_vect, frange(1));
hind = find_closest_ind(freq_vect, frange(2));

hars = [stimf(1):stimf(1):frange(2) stimf(2):-stimf(1):frange(1) stimf(2):stimf(1):frange(2)];
%       harmonics                   intermodulation              intermodulation pt 2
hars = unique(hars(hars>frange(1) & hars<frange(2)));

% need to check these, but it should be right
k = 1; % need to confirm this.
T = 2; 
halfrmband = (k+1)/(2*T); % remove the freq +- half the bandwidth

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