function subtightplotcleaner(fig, condfig, varargin)
% subtightplotcleaner: cleans up figure and axes in subtightplot plot
% 
% subtightplot(fig, condfig) takes the figure pointed to by the handle fig
% with the subplot with rows and colums pointed to by the two element
% vector condfig, and applies the default behavior (see below).
%
% subtightplot(fig, condfig, 'Name', 'Argument') takes name-argument pairs
% to control the default behaviour of subtightplotcleaner.
%
% The name-argument pairs and default behavious are described below.
%   heading:    string. Makes an overall heading on the figure using
%               suptitle if available. Alternately, download figtitle from
%               Matlab Central. Default is no heading.
%   catnames:   cell array of strings. names of the cats to be used for the
%               subplot headings. By default, this is empty. If set, it
%               labels the top row and first column.
%   cleanticks: logical. If set to true, cleans all subplot ticklabels. By
%               default, this is false, which leaves tickmarks on for the
%               bottom right-hand plot, moving the y-axis to the right.
%   topinds:    numeric vector. Defines the indicies used fron catnames to
%               label the subplots across the top row. Defaults to 18:34.
%   sideinds:   numeric vector. Same behaviour and defaults as topinds, but
%               acts on the subtitles along the first column. 
%   boxon:      logical. Default true. If false, the box (top and right
%               lines that delineate the plot) is turned off. Also sets
%               cleanticks to true.
%   xaxisscale: numeric or string. Sets behaviour of axis scaling. Use char
%               e.g. 'tight' to use a predefined behaviour or use a vector
%               to set limits manually. If using 'tight' or 'equal', you
%               must set this for both xaxisscale and yaxisscale.
%   yaxisscale: same as xaxisscale, but for the y axis.
%   xtickvals:  numeric vector. Sets XTick and XTickLabel property.

% Written by Rannee Li, Jul 2017
% Change log:
%   14 Jul 17:  Added description/help file. Added note about figtitle
%               function. 
%               Added axis stuff. 

% input parsing
p = inputParser;
addRequired(p, 'fig', @ishandle);

spccheck = @(x) isnumeric(x) && (numel(x)==2);
addRequired(p, 'condfig', spccheck)

charoremp = @(x) ischar(x) || isempty(x);
addParameter(p, 'heading', [], charoremp); % title over the whole thing

celloremp = @(x) iscell(x) || isempty(x);
addParameter(p, 'catnames', [], celloremp); % putting cat names on east and north axis

numorlog = @(x) isnumeric(x) || islogical(x);
addParameter(p, 'cleanticks', true, numorlog); % clear all ticks instead of only one
addParameter(p, 'topinds', 18:34, numorlog); % indicies for top names
addParameter(p, 'sideinds', 18:34, numorlog); % indices for side names
addParameter(p, 'box', true, numorlog);

scaleopts = @(x) isnumeric(x) || ischar(x);
addParameter(p, 'xaxisscale', [], scaleopts);
addParameter(p, 'yaxisscale', [], scaleopts);

cellornum = @(x) isnumeric(x) || iscell(x);
addParameter(p, 'xtickvals', [], cellornum)

parse(p, fig, condfig, varargin{:})

heading     = p.Results.heading;
catnames    = p.Results.catnames;
cleanticks  = p.Results.cleanticks;
topinds     = p.Results.topinds;
sideinds    = p.Results.sideinds;
boxon       = p.Results.box;
xaxisscale  = p.Results.xaxisscale;
yaxisscale  = p.Results.yaxisscale;
xtickvals   = p.Results.xtickvals;

% preclean setup
figure(fig)
nSubplot = prod(condfig);


% cycle through and perform functions

for n = 1:nSubplot
    subtightplot(condfig(1), condfig(2), n)
    
    % clean axis labels
    if ~boxon
        box(gca, 'off')
        % clean tick labels
        set(gca, 'XTickLabel', [], 'XTick', [])
        set(gca, 'YTickLabel', [], 'YTick', [])
    elseif ~cleanticks || (n ~= nSubplot && cleanticks)
        % clean tick labels
        if numel(xtickvals)>0
            set(gca, 'XTick', xtickvals)
        end
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
    else
        % only keep the bottom one on and move to RHS
        set(gca, 'YAxisLocation', 'right')
        if numel(xtickvals)>0
            set(gca, 'XTick', xtickvals)
            set(gca, 'XTickLabel', num2cell(xtickvals)')
        end
    end
    
    % add condition names
    % across the top
    if n <= condfig(2) && numel(catnames) > 0
        title(catnames{n}(topinds), 'interpret', 'none')
    end
    % on the side
    if mod(n, condfig(2))==1 && numel(catnames) > 0
        ylabel(catnames{n}(sideinds), 'interpret', 'none')
    end
    
    % scale if needed
    sptc_axscaler(gca, 'x', xaxisscale)
    sptc_axscaler(gca, 'y', yaxisscale)
    
end

% add overall heading if required
if numel(heading) > 0
    if exist('suptitle', 'file')==2
        % can call suptitle.
        suptitle(heading)
    else
        % call specialty function from matlab central
        figtitle(heading)
    end
end


function sptc_axscaler(ax, whichax, scaleropt)

if strcmpi(whichax, 'x')
    otherax = 'y';
elseif strcmpi(whichax, 'y')
    otherax = 'x';
else
    error('whichax value invalid')
end

if isempty(scaleropt)
    % empty actually breaks isnumeric, so just deal with it here
    return
elseif ischar(scaleropt)
    % get other axis values, set all, restore
    yl = get(ax, [otherax 'Lim']);
    axis(ax, scaleropt)
    set(ax, [otherax 'Lim'], yl)
elseif isnumeric(scaleropt)
    % set directly
    set(ax, [whichax 'Lim'], [scaleropt(1) scaleropt(end)])
    % this syntax avoids >2 element issue
end