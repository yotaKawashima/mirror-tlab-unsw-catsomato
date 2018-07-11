function data = call_biprref_v2(data_dir, filename, options, data)


if numel(data) < 1 % load data by calling the subfunction
    [datas, mvars, fname, singcondflag] = plothelp_loadopt(data_dir, filename, options);
    data.datas = datas;
    data.mvars = mvars;
    data.fname = fname;
    data.singcondflag = singcondflag;
end

%


figure(1); clf;
ch_label = draw_biprref_chlabfunction(data.mvars.label((prod(data.mvars.custom.spatialconfig)+1):end));

f_ind = find_closest_ind(data.mvars.freq{1}, options.frequency);
values = data.datas((prod(data.mvars.custom.spatialconfig)+1):end, f_ind, :);
if data.singcondflag
    % call single condition subfunction
    draw_biprref_sep(values, ch_label, data.mvars.custom.spatialconfig)
    colorbar()
    condname = options.condition{1};
else
    % call multicondition subfunction
    warning('multicondition plot not implemented yet')
    condname = 'allconditionsplot';
end

% -- print
for k = 1:numel(options.print)
    print(gcf, options.print{k}, [filename '_' condname '_' options.area '_bipolrrefplot'])
end



function [datas, mvars, fname, singcondflag] = plothelp_loadopt(data_dir, filename, options)
switch numel(options.condition)
    case 1
        % single condition
        fname = dir(fullfile(data_dir, [filename '*' options.area '*' options.condition{1} '*.mat'])); 
        singcondflag = true;
    case 0
        % all conditions
        fname = dir(fullfile(data_dir, [filename '*' options.condition '*.mat']));
        singcondflag = false;
    otherwise
        % invalid input
        error('Multiple conditions detected.')
end

numconds = numel(fname);

% -- loop contents pt 1
fprintf('Loading condition %2i of %2i\n', 1, numconds)
load(fullfile(data_dir, fname(1).name))

% -- processing for the first condition only
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
% -- loop contents, pt 2
avpow = mean(data.trial, 3);
datas(:, :, i) = avpow(:, fi1:fi2);

% actual loop
for i = 2:numconds
    fprintf('Loading condition %2i of %2i\n', i, numconds)
    load(fullfile(data_dir, fname(i).name))
    avpow = mean(data.trial, 3);
    datas(:, :, i) = avpow(:, fi1:fi2);
end

mvars.nconds = numconds;

