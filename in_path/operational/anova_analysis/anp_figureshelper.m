% helper file with templates for anp_figures.m

% Explanation of the fields
% 1. Select your type of plot
%   This is governed by the values in singfoi, singchan and classify. Each
%   of these is either true or false. Of the 8 options, only 5 are valid
%   and are listed below.
% 2. If needed, specify the frequencies and channels of interest
%   If singfoi is set to true, enter values for foi. Multiple frequencies
%   may be entered as a vector. If singchan is true when singfoi is false,
%   values in foi are used to draw vertical lines. For all other
%   circumstances, if singfoi is false then values in foi are ignored.
%   Similarly, if singchan is true, then use chan to specify the channels
%   of interest. If singchan is false, values in chan are ignored.
% 3. Set the p-value limits
%   If classify is true, then use classlim to specify the p-value for which
%   you want to set p<classlim as significant. If classify is false, any
%   value in this field is ignored.
%   If you want to use the p-value calculated using FDR, set classlim to
%   the string 'fdr', and modify the fdr field as required.  
%   For plots with colorbar, clim sets the maximum value of the colorbar.
%   Set clim to string 'fdr' if you want to use FDR to set the max value of
%   the colorbar, or use string 'default' if you want to use the full range
%   of values.
%   For plots with y-axis as p-value, elements in clim set the p-values for
%   which horizontal lines are drawn. If you would like to use FDR to find
%   1 or more of these values, set the FDR field. 
% 4. Set FDR
%   If you don't want to use FDR, set the fdr field to false or [].
%   Otherwise, fdr should be a 2x1 cell with the first element being the
%   Q-value and the second being the p type (either 'pID' or 'pN')
% TEMPLATES FOR DIFFERENT PLOTS FOLLOW
% remember that you can create an array of structs from an existing struct,
% e.g. b = a([1, 2]) is a 2x1 struct made from elements of a

% To plot p-values for all frequencies and all channels
options(1).singfoi      = false; % DO NOT EDIT
options(1).singchan     = false; % DO NOT EDIT
options(1).classify     = false; % DO NOT EDIT
% clim (below) is the max value of the colorbar. set to 'default' if you
% want to use default.
options(1).clim = 'fdr';
options(1).fdr = {q, 'pID'};


% To classify all channels and all frequencies                                
options(2).singfoi      = false; % DO NOT EDIT
options(2).singchan     = false; % DO NOT EDIT
options(2).classify     = true;  % DO NOT EDIT
% classlim (below) is the p-value for which p<classlim is significant.
options(2).classlim = 'fdr';
options(2).fdr = {q, 'pID'};

% To plot p-values for all frequencies for a single channel
options(3).singfoi      = false; % DO NOT EDIT
options(3).singchan     = true;  % DO NOT EDIT
options(3).classify     = false; % DO NOT EDIT
% clim (below) holds the values of horizontal dotted lines (p-values of
% interest). If you don't want them, set to [].
options(3).clim = [0.05, 0.01];
options(3).fdr = {q, 'pID'}; % set this to 'false' if you want to clim only
% foi (below) holds the values of the vertical dotted lines (frequencies of
% interest). If you don't want them, set to [].
options(3).foi = [23, 46, 200];
% chan (below) holds the channels. Each channel appears on a separate plot.
options(3).chan = [178, 146];

% To plot p-values for all channels at a given frequency
options(4).singfoi      = true;  % DO NOT EDIT
options(4).singchan     = false; % DO NOT EDIT
options(4).classify     = false; % DO NOT EDIT
% clim (below) is the max value of the colorbar. set to [] if you want to
% use default.
options(4).clim = 'fdr';
options(4).fdr = {q, 'pID'}; 
% foi (below) are the frequencies of interest. Each frequency has its own
% plot.
options(4).foi = [23, 200];

% To classify channels at a given frequency
options(5).singfoi      = true;  % DO NOT EDIT
options(5).singchan     = false; % DO NOT EDIT
options(5).classify     = true;  % DO NOT EDIT
% foi (below) are the frequencies of interest. Each frequency has its own
% plot.
options(5).foi = [23, 200];
% classlim (below) is the p-value for which p<classlim is significant.
options(5).classlim = 'fdr';
options(5).fdr = {q, 'pID'};