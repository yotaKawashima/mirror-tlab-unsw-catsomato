function data_out = highgammapower_y(data_dir, loadname, fband, foi, hbw)

% highgammapower: calculate the total power in the high gamma band
% 
%   data_out = highgammapower(data_dir, loadname, fband, foi, hbw) loads
%   the data from data_dir with the names indicated by loadname. fband
%   describes the band of frequencies, with foi being the frequencies to be
%   removed (response frequencies, harmonics and intermodulation). hbw is
%   the half-bandwith to be removed around foi.
%
%   data.trial becomes a nCond-1 x 1 cell array, with each cell containing
%   a nChan x nFreq x nTrial matrix. 

% Written by Yota, April 2020 
% Updata Rennee's code. 
% Just extract VELogP data in High Gamma Band. Foi (i.e. tagged
% frequencies) are not included. At this point, mean across frequnices and
% trials is not taken.

% load the first dataset
% start the 1. In Rennee's code, number 1 is removed.
load(fullfile(data_dir, loadname(1).name))

% make mask
fs = data.freq{1}';

[mask, maskedf] = hgp_makemask(fs, foi, fband, hbw);

data_out = data;
data_out.freq{1} = maskedf;
t_out = cell(numel(loadname), 1);
data_out.trial = t_out;

% find and label data
data_out.trial{1} = data.trial(:, logical(mask), :);
data_out.datalabels{1} = loadname(1).name(18:34);

% do the rest of the names
for k = 2:numel(loadname)
    
    % load
    load(fullfile(data_dir, loadname(k).name))
    
    % find and label data
    data_out.trial{k} = data.trial(:, logical(mask), :);
    data_out.datalabels{k} = loadname(k).name(18:34);

end


% update filename
data_out.custom.filename = [data.custom.filename '_hgpcomp'];
data_out.custom.filename(18:34) = 'FxxxAxxx_FxxxAxxx';


% ----- SUBFUNCTIONS FOLLOW -----
function [mask, maskedf] = hgp_makemask(fs, foi, fband, hbw)
% create problem band mask

mask = ones(size(fs));

% mask upper and lower
[~, m] = find_closest(fs, fband(1));
mask(1:m-1) = 0;

[~, m] = find_closest(fs, fband(2));
mask(m+1:end) = 0;

% how many samples in a halfbandwidth?
spacing = mean(diff(fs));
hbwS = round(hbw/spacing);

% mask foi +- halfbandwidth
for k = 1:numel(foi)
    [~, m] = find_closest(fs, foi(k));
    
    mask(m-hbwS:m+hbwS) = 0;
end

maskedf = fs(logical(mask));




