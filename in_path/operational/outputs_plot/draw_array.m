function H = draw_array(data, cfg)

% draws data in array format
% called by call_plots.m
% returns handle H if nargout > 0

% find only the freq of interest
[~, ind] = find_closest(data.freq{1}, cfg.foi);

% preallocate the shape
nNonrref = data.custom.spatialconfig(1) * data.custom.spatialconfig(2);
config = zeros(data.custom.spatialconfig);

% fill the shape
config(:) = data.trial(1:nNonrref, ind);

% draw
imagesc([], [], real(config), [cfg.limits(1) cfg.limits(2)])

axis ij

if nargout < 1
    clear H
end