function [mask, maskedf] = makefreqmask(fs, foi, fband, hbw)

% makefreqmask: makes a mask to remove tagged frequencies.
%
% mask = makefreqmask(fs, foi) returns a logical the same length as fs
% where values equal to foi are equal to zero. All other vallues are 1.
%
% mask = makefreqmask(fs, foi, fband) also sets values less than fband(1)
% and greater than fband(end) equal to 0.
%
% mask = makefreqmask(fs, foi, fband, hbw) also sets values within hbw of
% the fois to 0.
%
% [mask, maskedf] = makefreqmask(...) returns the masked values of fs

% set defaults
if nargin < 3
    fband(1) = fs(1);
    fband(2) = fs(end);
end
if nargin < 4
    hbw = 0;
end

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
    
    mask(max([m-hbwS,0]):min([m+hbwS, numel(fs)])) = 0;
end

mask = logical(mask);

if nargout > 1
    maskedf = fs(mask);
end