function [mask, maskedf] = makefreqmask(fs, foi, fband, hbw)

% makefreqmask: makes a mask to remove tagged frequencies.
%
% mask = makefreqmask(fs, foi) returns a vector the same length as fs where
% values equal to foi are equal to zero. All other vallues are 1.
%
% mask = makefreqmask(fs, foi, fband) also sets values less than fband(1)
% and greater than fband(end) equal to 0.
%
% mask = makefreqmask(fs, foi, fband, hbw) also sets values within hbw of
% the fois to 0.
%
% [mask, maskedf] = makefreqmask(...) returns the masked values of fs



% create problem band mask

mask = ones(size(fs));

% mask upper and lower
[~, m] = find_closest(fs, fband(1));
mask(1:m-1) = 0;

[~, m] = find_closest(fs, fband(end));
mask(m+1:end) = 0;

% how many samples in a halfbandwidth?
spacing = mean(diff(fs));
hbwS = floor(hbw/spacing)+1;

% mask foi +- halfbandwidth
for k = 1:numel(foi)
    [~, m] = find_closest(fs, foi(k));
    
    mask(m-hbwS:m+hbwS) = 0;
end



if nargout > 1
    maskedf = fs(logical(mask));
end