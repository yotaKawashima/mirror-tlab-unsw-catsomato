% create gatherdata cfg

% if all conds:
nConds = 16;

CatNum   = 'C20110808';
RecNum   = 'R03';
StimType = 'TStim';
Area     = 'S1';
DataType = 'epoched_rsampsl_biprref_evkresp_pwrspec_f0to250_meantrl'; %_newzscr_zscored_dwnsmpl

CondNames = {'F023A000_F200A000', 'F023A000_F200A004', 'F023A000_F200A007', 'F023A000_F200A016', ...
             'F023A040_F200A000', 'F023A040_F200A004', 'F023A040_F200A007', 'F023A040_F200A016', ...
             'F023A079_F200A000', 'F023A079_F200A004', 'F023A079_F200A007', 'F023A079_F200A016', ...
             'F023A159_F200A000', 'F023A159_F200A004', 'F023A159_F200A007', 'F023A159_F200A016'};

cfg.filename = cell(1, nConds);
for i = 1:nConds
    cfg.filename{i} = [CatNum '_' RecNum '_' StimType '_' Area '_' ...
        CondNames{i} '_' DataType '.mat'];
end



cfg.datadir = ['/Users/ranneeli/Documents/MATLAB/ftformat_pipe/data/' DataType '/'];



% set channels and trials
% leave empty to keep all.
cfg.channels = [];
cfg.trials   = [];

% cfg.samples  = [];
% 
% foi = 23;
% if ~isempty(foi)
%     [~, cfg.samples] = findclosest(;
% end

drawnow;

% %%
% % fprintf('[mixed, plotcfg] = gatherdata(cfg);\n')
% tic;
[mixed, plotcfg] = gatherdata(cfg);
% save -v7.3 testing mixed plotcfg
% time_take = toc %#ok<NOPTS>
% fprintf('Mixed\n')
% 
% %%
% 
% mixedmean = mixed;
% for i = 1:nConds
%     fprintf('Condition %02i\n', i)
%     tic
%     
%     mixedmean(i).trial = mean(mixedmean(i).trial, 3);
%     
%     taken = toc;
%     fprintf('\tTime = %4.2f sec\n', taken)
% end

if strcmp(Area, 'S1')
    plotcfg.biprref_index = 101:280; %65:176; %
elseif strcmp(Area, 'S2')
    plotcfg.biprref_index = 65:176; %101:280; %
else
    error('Plz.')
end
% plotcfg.foi           = 200;

%%
