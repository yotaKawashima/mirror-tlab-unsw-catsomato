function data_out = bipolreref(data)

% Bipolar rereferencing in both horizontal and vertical directions.
%
%    This function keeps the original data as well and creates a new set of
%    channel labels.
%    Data is assumed to have been collected in a rectangular array format.

% Written by Rannee Li, Dec 2014

% copy everything
data_out = data;

% extract the following data parameters
num_rawdat    = data.custom.nsignals;
ConfigRows    = data.custom.spatialconfig(1);
ConfigColumns = data.custom.spatialconfig(2);
nSamples      = data.custom.nsamples;
nTrials       = data.custom.ntrials;

% calculate conditions
num_rerefv = (ConfigColumns-1) * ConfigRows;
num_rerefh = (ConfigRows-1) * ConfigColumns;
totalnames = num_rawdat + num_rerefv + num_rerefh;

% preallocate outputs
chlabels = cell(totalnames, 1);
biprref = zeros(totalnames, nSamples, nTrials);

% rename matrix
matrix = data.trial;

% remove NaNs
matrix(isnan(matrix)) = 0;

% define arrayconfig
arrayconfig = zeros(ConfigRows, ConfigColumns);
arrayconfig(:) = 1:num_rawdat;

%% raw data

for i = 1:num_rawdat
    % name
    chlabels{i} = data_out.label{i};
    
end

% copy data
biprref(1:num_rawdat, :, :) = matrix;

%% vertical rereference

for confc = 1:ConfigRows
    % find numbers of all the channels in the column
    col_inds = arrayconfig(:, confc);
    
    % find the channel numbers for the subtraction
    sub_ind1 = col_inds(2:end);
    sub_ind2 = col_inds(1:end-1);
    
    % collect the values that are going to subtract
    sub_arg1 = matrix(sub_ind1, :, :);
    sub_arg2 = matrix(sub_ind2, :, :);
    
    % maths
    refd = sub_arg1 - sub_arg2;
    
    % move to reref
    ref_ind = num_rawdat + (confc-1)*(ConfigRows-1) + 1;
    biprref(ref_ind:(ref_ind + ConfigRows -2), :, :) = refd;
    
    chnum = 1;
    for i = ref_ind:(ref_ind + ConfigRows - 2)
        chlabels{i} = ['rerefv_ch' num2str(sub_ind1(chnum), '%03i') '-ch' num2str(sub_ind2(chnum), '%03i')];
        
        % increment counter
        chnum = chnum + 1;
    end
    
end

%% horizontal rereference

for confr = 1:ConfigColumns
    % find the channel numbers of all channels in the row
    row_inds = arrayconfig(confr, :);
    
    % find the channel numbers for the subtraction
    sub_ind1 = row_inds(2:end);
    sub_ind2 = row_inds(1:end-1);
    
    % collect the values that are going to subtract
    sub_arg1 = matrix(sub_ind1, :, :);
    sub_arg2 = matrix(sub_ind2, :, :);
    
    % maths
    refd = sub_arg1 - sub_arg2;
    
    % move to reref
    ref_ind = num_rawdat + num_rerefv + (confr-1)*(ConfigColumns-1) + 1;
    biprref(ref_ind:(ref_ind + ConfigColumns -2), :, :) = refd;
    
    chnum = 1;
    for i = ref_ind:(ref_ind + ConfigColumns - 2)
        chlabels{i} = ['rerefh_ch' num2str(sub_ind1(chnum), '%03i') '-ch' num2str(sub_ind2(chnum), '%03i')];
        
        % increment counter
        chnum = chnum + 1;
    end
    
end

%% struc allocation
data_out.trial = biprref;
data_out.label = chlabels;
data_out.custom.nsignals = totalnames;

data_out.custom.filename = [data_out.custom.filename '_biprref'];

data_out.custom.datatype.isbiprref = true;