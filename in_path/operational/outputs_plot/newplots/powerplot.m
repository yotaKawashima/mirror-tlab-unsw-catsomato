function [data, axhand, fighand] = powerplot(data_dir, filename, options, data)

% plots the power across all of the conditions
% 
% Usage:
%   data = powerplot(data_dir, filename, options, data) plots the stimulus 
%
%   [data, axhand, fighand] = powerplot(...) returns the axis handle in
%   axhand and the figure handle in fighand.
%
%
% Fields in options:
%   frange: frequency range of interest [fstart, fend]. If there are
%       more than two elements, the first and last are used.
%   area: the area that the data is coming from. 'S1' or 'S2'
%   channel: the channel number. double int.
%   version: the matlab version (for setting axis titles). The year as an
%       integer.
%   print: cell array of strings used in the print statement as the file
%   	type. 


if numel(data)<1
    % load data
    [datas, mvars, fname] = pp_loaddata(data_dir, filename, options);
    data.datas = datas;
    data.mvars = mvars;
    data.fname = fname;
end

ylims = find_lims(data.datas);

set(0,'DefaultTextInterpreter','none');



x2 = data.mvars.custom.subplotconfig(2);
x1 = data.mvars.custom.subplotconfig(1);

figure(2);clf
for k = 1:data.mvars.nconds
    subtightplot(x2, x1, k)
    plot(data.mvars.freq{1}, data.datas(options.channel, :, k))
    gca_set('XLim', options.frange, options.version)
    gca_set('YLim', ylims, options.version)
    
    if k <= x2 % top row
        [t, ~] = regexp(data.fname(k).name, 'F...A...', 'match', 'split');
        title(t{2}, 'FontWeight', 'Normal')
    end
    
    if mod(k, x1) == 0 % right column
        [t, ~] = regexp(data.fname(k).name, 'F...A...', 'match', 'split');
        ylabel(t{1})
        gca_set('yaxislocation','right', 2010); % having issues???
    end
    
    if k ~= x2*(x1-1)+1 % not bottom left
        gca_set('XTickLabel', [], options.version)
        gca_set('YTickLabel', [], options.version)
        
    else % bottom left
        xlabel('frequency (Hz)')
        ylabel('log power (dB)')
    end
end

% print
for k = 1:numel(options.print)
    print(gcf, options.print{k}, [filename '_' options.area '_CH' ...
        num2str(options.channel, '%03i') '_powerspecplot'])
end

if nargout > 1
    fighand = gcf;
    axhand = gca;
end



function [datas, mvars, fname] = pp_loaddata(data_dir, filename, options)

fname = dir(fullfile(data_dir, [filename, '*' options.area '*.mat'])); % always working in S2

numconds = numel(fname);

fprintf('Loading condition %2i of %2i\n', 1, numconds)
load(fullfile(data_dir, fname(1).name))

% some processing
if isfield(options, 'frange')
    fs = sort([options.frange(1) options.frange(end)]);
    fi1 = find_closest_ind(data.freq{1}, fs(1));
    fi2 = find_closest_ind(data.freq{1}, fs(2)); 
    
else
    fi1 = 1; % first frequency index
    fi2 = size(data.trial, 2);
end

% place into variables
mvars = rmfield(data, 'trial');
mvars.freq{1} = mvars.freq{1}(fi1:fi2);
[a, ~, ~] = size(data.trial);
datas = zeros(a, fi2-fi1+1, numconds);

i = 1;
% loop contents
avpow = mean(data.trial, 3);
datas(:, :, i) = avpow(:, fi1:fi2);

for i = 2:numconds
    fprintf('Loading condition %2i of %2i\n', i, numconds)
    load(fullfile(data_dir, fname(i).name))
    avpow = mean(data.trial, 3);
    datas(:, :, i) = avpow(:, fi1:fi2);
end

mvars.nconds = numconds;