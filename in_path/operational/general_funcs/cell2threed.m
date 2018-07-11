function data_out = cell2threed(data_in, nChannels, nSamples, nTrials)

% Converts data_in from a cell array to a 3d matrix data_out
% Written by Rannee Li, Jan 2015
%
% dimensions of data_out:    nChannels x nSamples x nTrials
%
% dimensions of data_in:   1 x nTrials cell, with each cell a double with
%                           dimensions nChannels x nSamples

data_out = zeros(nChannels, nSamples, nTrials);

for i = 1:nTrials
    data_out(:, :, i) = data_in{i};
end