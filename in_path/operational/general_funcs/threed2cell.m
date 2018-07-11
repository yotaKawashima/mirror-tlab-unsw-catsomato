function data_out = threed2cell(data_in, nTrials)

% Converts data_in from a 3d matrix to a cell array data_out
% Written by Rannee Li, Jan 2015
%
% dimensions of data_in:    nChannels x nSamples x nTrials
%
% dimensions of data_out:   1 x nTrials cell, with each cell a double with
%                           dimensions nChannels x nSamples

data_out = cell(1, nTrials);

for i = 1:nTrials
    data_out{i} = squeeze(data_in(:, :, i));
end