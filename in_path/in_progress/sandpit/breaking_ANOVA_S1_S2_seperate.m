% Plot ANOVA results(_datout) cross all experiments with S1 and S2
% seperately, and see meeting note 02011300 1,2-1
% Wrriten by Nao, editted by Phoebe, 01 Feb 2016

% This code loads ANOVA results('_datout') for all experiments, concatenate
% all channels (S1, S2 seperately), plot p values analysis with
% frequencies range 0-220
clear;
clc;
close all;
dbstop if error

Plot.Concatenate = 0;
Plot.Mean = 0;
Plot.Median = 0;
Plot.Prctile = 0;

ReponseFreq = 1;
%% Data need to be editted
% Pick S1 or S2?
A = {'_S1'};
% A = {'_S2'};
switch A{:}
    case {'_S1'}
        nAllChanPerExp = 180;
        UnipoChan = 100;
    case {'_S2'}
        nAllChanPerExp = 112;
        UnipoChan = 64;
end
TargFreq = 23;  % Pick 23Hz or 200 Hz,(as the response frequency)
% TargFreq = 200;
%% Prepare data for all experiments (note: pay attention to the total number of expperiments)
filename = ['C20110511_R02_TStim';'C20110510_R05_TStim';...
    'C20110510_R06_TStim';'C20110808_R01_TStim';'C20110808_R03_TStim';...
    'C20110808_R04_TStim';'C20110808_R06_TStim';'C20110808_R09_TStim';'C20110808_Rx4_TStim'];
% disp(size(filename));
func_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/';
addpath(genpath(func_dir))
data_dir = '/media/phoebeyou/My Passport/Spencers_Cat_Data/data/epoched_rsampsl_biprref_evkresp_cmtspwr_adatout/';

% def_name = {' 23 Hz dependence', ' 200Hz dependence', 'interaction'};

figure(1)
set(gcf,'visible','off')

all_pval =[];
for iExp = 1:size(filename,1)
    for iRecordSite = 1:numel(A)
        %    load anova results
        %         disp(filename(iExp,:))
        filename_header = [filename(iExp,:) A{iRecordSite}];
        fname = dir(fullfile(data_dir,[filename_header '*_P*_anov*_adatout.mat']));
        load([data_dir fname.name],'pvals','metavars')
        %         test_pval(iExp,:) = pvals;
        %         load([data_dir fname.name])
        anovaResult{iExp,iRecordSite} = metavars; %data;
        all_pval = [all_pval permute(pvals,[2 1 3])];
        %         disp(size(all_pval))
        %         find(isnan(all_pval(1,:,1))')
        if any(anovaResult{iExp,iRecordSite}.freq{1} ~= anovaResult{1,1}.freq{1} )
            keyboard
        end
        
    end
end
% Remove the NaN channels
goodChan = find(~isnan(all_pval(1,:,1))');
good_pval = all_pval(:,goodChan,:);
% size(all_pval)
% size(good_pval)
%% ***************************************************************************************************************************
%% Plot all expriments (after concatenate S1&S2 individually)
% depedence on 23Hz, 200Hz and interaction
if Plot.Concatenate
    figure(1),clf
    set(gcf,'visible','on')
    imagesc(anovaResult{iExp,iRecordSite}.freq{1},[], log10(good_pval(:, :, 1))',[-5 0])
    h = colorbar;
    title(h,'log10(p value)')
    title(['p values in dependence of 200Hz about ' A],'interpret','none')
    xlabel('Frequency(Hz)')
    ylabel('Channel')
    print(1, '-dpng', 'pval_depend200' )
    
    figure(2),clf
    set(gcf,'visible','on')
    imagesc(anovaResult{iExp,iRecordSite}.freq{1},[], log10(good_pval(:, :, 2))',[-5 0])
    h = colorbar;
    title(h,'log10(p value)')
    title(['p values in dependence of 23Hz about ' A],'interpret','none')
    xlabel('Frequency(Hz)')
    ylabel('Channel')
    print(2, '-dpng', 'pval_depend23_')
    
    figure(3),clf
    set(gcf,'visible','on')
    imagesc(anovaResult{iExp,iRecordSite}.freq{1},[], log10(good_pval(:, :, 3))',[-5 0])
    h = colorbar;
    title(h,'log10(p value)')
    title(['p values in dependence of interaction about ' A],'interpret','none')
    xlabel('Frequency(Hz)')
    ylabel('Channel')
    print(3, '-dpng', 'pval_dependInteraction_' )
end
% Plot mean
if Plot.Mean
    figure(4),clf
    hold on
    plot(anovaResult{iExp,iRecordSite}.freq{1},mean(log10(good_pval(:, :, 1)')),'r') % 200
    plot(anovaResult{iExp,iRecordSite}.freq{1},mean(log10(good_pval(:, :, 2)')),'g') % 23
    plot(anovaResult{iExp,iRecordSite}.freq{1},mean(log10(good_pval(:, :, 3)')),'b') % interaction
    title(['Mean for all experiments about ' A],'interpret','none')
    ylabel('log10(p value)')
    xlabel('Frequency(Hz)')
    legend('200Hz','23Hz','Interaction','Location','SouthOutside','Orientation','horizontal')
    print(4, '-dpng', 'Mean_AllExp_' )
end
% Plot median
if Plot.Median
    figure(5),clf
    hold on
    plot(anovaResult{iExp,iRecordSite}.freq{1},median(log10(good_pval(:, :, 1)')),'r') % 200
    plot(anovaResult{iExp,iRecordSite}.freq{1},median(log10(good_pval(:, :, 2)')),'g') % 23
    plot(anovaResult{iExp,iRecordSite}.freq{1},median(log10(good_pval(:, :, 3)')),'b') % interaction
    title(['Median for all experiments about ' A],'interpret','none')
    ylabel('log10(p value)')
    xlabel('Frequency(Hz)')
    legend('200Hz','23Hz','Interaction','Location','SouthOutside','Orientation','horizontal')
    print(5, '-dpng', 'Median_AllExp_' )
end
% Plot percentiles with a percentage (interested in percentage of 0.01, 0.02, 0.05, 0.5)
if Plot.Prctile
    percen = 5/100;
    figure(6),clf
    hold on
    plot(anovaResult{iExp,iRecordSite}.freq{1},prctile(log10(good_pval(:, :, 1)'),percen),'r') % 200
    plot(anovaResult{iExp,iRecordSite}.freq{1},prctile(log10(good_pval(:, :, 2)'),percen),'g') % 23
    plot(anovaResult{iExp,iRecordSite}.freq{1},prctile(log10(good_pval(:, :, 3)'),percen),'b') % interaction
    title(['Percentile for all experiments percentage of 0.05, about ' A],'interpret','none')
    ylabel('log10(p value)')
    xlabel('Frequency(Hz)')
    legend('200Hz','23Hz','Interaction','Location','SouthOutside','Orientation','horizontal')
    % print(6,'-dpng','Percentile0.05_AllExp_S2' )
end

%% 2-1 ****************************************************************************************************************************
%% find channels that are significnat for main effect of 23Hz stim only
% not for 200Hz stim and interaction
freq = anovaResult{iExp,iRecordSite}.freq{1};

%% at 23Hz response

if ReponseFreq
    fig7 = 0;
    fig8 = 0;
    fig9 = 1;
    figPie = 0;
    %% The effect factor a for 23, b for 200, c for interaction (as the stimulus)
    Effect = ['  a';'  b';'  c';' ab';' ac';' bc';'abc'];
    
    [~,iTargFreq]= min(abs(freq - TargFreq));
    disp([num2str(freq(iTargFreq)) ' Hz response'])
    % Create a matrix to contain the number of channles satisfy
    % Effect(iEffect), and then plot the distribution
    ChanDistribution = zeros(1,length(Effect));
    for iEffect =  1:length(Effect)
        switch Effect(iEffect,:)
            case {'  a'}
                vAllChan = find(all_pval(iTargFreq,:,2) <0.05 & all_pval(iTargFreq,:,1) > 0.05 & all_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <a> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,2);
            case {'  b'}
                vAllChan = find(all_pval(iTargFreq,:,2) >0.05 & all_pval(iTargFreq,:,1) < 0.05 & all_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <b> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,1);
            case {'  c'}
                vAllChan = find(all_pval(iTargFreq,:,2) >0.05 & all_pval(iTargFreq,:,1) > 0.05 & all_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <c> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,3);
            case {' ab'}
                vAllChan = find(all_pval(iTargFreq,:,2) <0.05 & all_pval(iTargFreq,:,1) < 0.05 & all_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <ab> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,1) .* all_pval(iTargFreq,vAllChan,2);
            case {' ac'}
                vAllChan = find(all_pval(iTargFreq,:,2) <0.05 & all_pval(iTargFreq,:,1) > 0.05 & all_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <ac> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,2) .* all_pval(iTargFreq,vAllChan,3);
            case {' bc'}
                vAllChan = find(all_pval(iTargFreq,:,2) >0.05 & all_pval(iTargFreq,:,1) < 0.05 & all_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <bc> effect only is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,1) .* all_pval(iTargFreq,vAllChan,3);
            otherwise
                vAllChan = find(all_pval(iTargFreq,:,2) <0.05 & all_pval(iTargFreq,:,1) < 0.05 & all_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <abc> effect is ' num2str(length(vAllChan))])
                pval_targFreq = all_pval(iTargFreq,vAllChan,1) .* all_pval(iTargFreq,vAllChan,2) .* all_pval(iTargFreq,vAllChan,3);
        end
        ChanDistribution(1,iEffect) = length(vAllChan);
        %         switch Effect(iEffect,:)
        %             case {'  a'}
        %                 vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) > 0.05);
        %                 disp(['The number of electrons satisfy <a> effect only is ' num2str(length(vChan))])
        %             case {'  b'}
        %                 vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) > 0.05);
        %                 disp(['The number of electrons satisfy <b> effect only is ' num2str(length(vChan))])
        %             case {'  c'}
        %                 vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) < 0.05);
        %                 disp(['The number of electrons satisfy <c> effect only is ' num2str(length(vChan))])
        %             case {' ab'}
        %                 vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) > 0.05);
        %                 disp(['The number of electrons satisfy <ab> effect only is ' num2str(length(vChan))])
        %             case {' ac'}
        %                 vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) < 0.05);
        %                 disp(['The number of electrons satisfy <ac> effect only is ' num2str(length(vChan))])
        %             case {' bc'}
        %                 vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) < 0.05);
        %                 disp(['The number of electrons satisfy <bc> effect only is ' num2str(length(vChan))])
        %             otherwise
        %                 vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) < 0.05);
        %                 disp(['The number of electrons satisfy <abc> effect is ' num2str(length(vChan))])
        %         end
        
    end
    for iEffect =  1 %1:length(Effect)       
        
        %% Find minimum pval and the location of it in the original concatenate
        % matrix(all_pval), further determine which experiment and channel it is
        [min_pval,iMin_pval] = min(pval_targFreq);
        iMin_pval_chan = vAllChan(iMin_pval);
        
        %        nAllChanPerExp = nAllChanS1PerExp + nAllChanS2PerExp;
        %        iRealChanPerExp = mod(iMin_pval_chan - 1 ,nAllChanPerExp ) + 1
        iRealChanPerExp = mod(iMin_pval_chan - 1 ,nAllChanPerExp ) + 1;
        
        iRealExp = ceil ( iMin_pval_chan / nAllChanPerExp );
        disp('The minimum pval & channel after concatenating is ')
        disp([num2str(min_pval) ' : ' num2str(iMin_pval_chan)])
        disp(['in Exp: ' num2str(filename(iRealExp,:)) ' : chan = ' num2str(iRealChanPerExp)])
        
        %% After find the matched exp and channel, load all conditions for this exp at this channel
        MatchFile = [filename(iRealExp,:) A{:}];
        dir_sec = [data_dir(1:end-48) 'epoched_rsampsl_biprref_evkresp_cmtspwr/'];
        loadname = dir(fullfile(dir_sec, [MatchFile '*.mat']));
        % Define num of conditions for this exp
        switch filename(iRealExp,:)
            case {'C20100414_R09_TStim','C20110511_R02_TStim'}
                ConValue23 = ['000' ;'019'; '040'; '079'; '159'];
                ConValue200 = ['000'; '002'; '004' ;'007';'015'];
            case {'C20110510_R05_TStim','C20110510_R06_TStim'}
                ConValue23 = ['000'; '010'; '020'; '040'; '079'; '159'];
                ConValue200 = ['000'; '001'; '002'; '004'; '007'; '015'];
            case {'C20110808_R01_TStim','C20110808_R03_TStim','C20110808_R04_TStim'}
                ConValue23 = ['000'; '040'; '079'; '159'];
                ConValue200 = ['000'; '004'; '007'; '016'];
            otherwise % {'C20110808_R06_TStim','C20110808_R09_TStim','C20110808_Rx4_TStim'}
                ConValue23 = ['000'; '019'; '040'; '079'; '159'];
                ConValue200 = ['000'; '004'; '007'; '016'; '031'];
        end
        xlength = length(ConValue200);
        ylength = length(ConValue23);
        %% Load power spectrum result at channel no. iRealChanPerExp+UnipoChan
        %         load(dir_sec, [(filename(iRealExp,1:end-5)) 'Chan' iRealChanPerExp]);
        %                 for iFile = 1:length(loadname)
        %         %                     load(fullfile(dir_sec, loadname(iFile).name));
        %                             disp(iFile)
        %                             disp(loadname(iFile).name)
        %                             if iFile == 1
        %                                 allData = zeros(length(loadname),size(data.trial,2),size(data.trial,3));
        %                             end
        %                             % skip first 100 or 64 unipolar
        %                             allData(iFile,:,:)= data.trial(iRealChanPerExp+UnipoChan,:,:);  %100 for S1, 64 for S2
        %                         end
        
        %% Load anova result to check
        %  [a,b] = min(all_pval(iTargFreq,vAllChan,2));
        switch Effect(iEffect,:)
            case {'  a'}
                vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <a> effect only is ' num2str(length(vChan))])
                %                 [a,b] = min(good_pval(iTargFreq,vChan,2));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,2));
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,2));
                %                 pval_allChan = all_pval(iTargFreq,:,2);
                pval_allChan = good_pval(iTargFreq,:,2);
            case {'  b'}
                vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <b> effect only is ' num2str(length(vChan))])
                %                 [a,b] = min(good_pval(iTargFreq,vChan,1));
                %                 [minPval_qualifying,iMinPval_qualifying]  = min(all_pval(iTargFreq,vAllChan,1));
                %                 pval_allChan = all_pval(iTargFreq,:,1);
                [minPval_qualifying,iMinPval_qualifying]  = min(good_pval(iTargFreq,vChan,1));
                pval_allChan = good_pval(iTargFreq,:,1);
            case {'  c'}
                vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <c> effect only is ' num2str(length(vChan))])
                %                 [a,b] = min(good_pval(iTargFreq,vChan,3));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,3));
                %                 pval_allChan = all_pval(iTargFreq,:,3);
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,3));
                pval_allChan = good_pval(iTargFreq,:,3);
            case {' ab'}
                vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) > 0.05);
                disp(['The number of electrons satisfy <ab> effect only is ' num2str(length(vChan))])
                %                 [a,b] = min(good_pval(iTargFreq,vChan,1).*good_pval(iTargFreq,vChan,2));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,2).*all_pval(iTargFreq,vAllChan,1));
                %                 pval_allChan = all_pval(iTargFreq,:,1) .* all_pval(iTargFreq,:,2);
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,2).*good_pval(iTargFreq,vChan,1));
                pval_allChan = good_pval(iTargFreq,:,1) .* good_pval(iTargFreq,:,2);
            case {' ac'}
                vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) > 0.05 & good_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <ac> effect only is ' num2str(length(vChan))])
                %                 [a,b] = min(good_pval(iTargFreq,vChan,2).*good_pval(iTargFreq,vChan,3));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,2).*all_pval(iTargFreq,vAllChan,3));
                %                 pval_allChan = all_pval(iTargFreq,:,2) .* all_pval(iTargFreq,:,3);
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,2).*good_pval(iTargFreq,vChan,3));
                pval_allChan = good_pval(iTargFreq,:,2) .* good_pval(iTargFreq,:,3);
            case {' bc'}
                vChan = find(good_pval(iTargFreq,:,2) >0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <bc> effect only is ' num2str(length(vChan))])
                %                 [a,b]= min(good_pval(iTargFreq,vChan,1).*good_pval(iTargFreq,vChan,3));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,1).*all_pval(iTargFreq,vAllChan,3));
                %                 pval_allChan = all_pval(iTargFreq,:,1) .* all_pval(iTargFreq,:,3);
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,1).*good_pval(iTargFreq,vChan,3));
                pval_allChan = good_pval(iTargFreq,:,1) .* good_pval(iTargFreq,:,3);
            otherwise
                vChan = find(good_pval(iTargFreq,:,2) <0.05 & good_pval(iTargFreq,:,1) < 0.05 & good_pval(iTargFreq,:,3) < 0.05);
                disp(['The number of electrons satisfy <abc> effect is ' num2str(length(vChan))])
                %                 [a,b]= min(good_pval(iTargFreq,vChan,1).*good_pval(iTargFreq,vChan,2).*good_pval(iTargFreq,vChan,3));
                %                 [minPval_qualifying,iMinPval_qualifying] = min(all_pval(iTargFreq,vAllChan,1).*all_pval(iTargFreq,vAllChan,2).*all_pval(iTargFreq,vAllChan,3));
                %                 pval_allChan = all_pval(iTargFreq,:,1) .* all_pval(iTargFreq,:,2) .* all_pval(iTargFreq,:,3);
                [minPval_qualifying,iMinPval_qualifying] = min(good_pval(iTargFreq,vChan,1).*good_pval(iTargFreq,vChan,2).*good_pval(iTargFreq,vChan,3));
                pval_allChan = good_pval(iTargFreq,:,1) .* good_pval(iTargFreq,:,2) .* good_pval(iTargFreq,:,3);
        end
        if min_pval == minPval_qualifying
            disp('Match.')
        else
            disp('Calculated min pval mismatch with real min p val.');
        end
        
        disp('-------------------------------------------------------------------------------------')
        %% Plot main effect of Effect(iEffect),how 23,200,interaction looks like at 23Hz reponse
        if fig7
            figure(7),clf
            %             plot(log10(good_pval(iTargFreq,vChan,2)),'g'); % 23
            plot(log10(good_pval(iTargFreq,vChan,2)),'g'); % 23
            hold on
            %             plot(log10(good_pval(iTargFreq,vChan,1)),'r'); % 200
            plot(log10(good_pval(iTargFreq,vChan,1)),'r'); % 200
            %             plot(log10(good_pval(iTargFreq,vChan,3)),'b'); % interaction
            plot(log10(good_pval(iTargFreq,vChan,3)),'b'); % interaction
            plot(iMinPval_qualifying*[1 1], ylim,'k--')
            ylabel('Significant log(pval)','FontSize',12)
            xlabel('Channel','FontSize',12)
            legend('23','200','interaction',['minPval:' num2str(iMinPval_qualifying)],'Location','SouthOutside','Orientation','Horizontal')
            title( ['S1 ' Effect(iEffect,:) ' main effect at ' num2str(TargFreq) 'Hz response,NumChan:' num2str(length(vAllChan)) ],'Fontsize',15);
            print(7, '-dpng', [Effect(iEffect,:) ' main efect response'])
        end
        %%
        if fig8
            figure(8),clf
            hold on
            plot(log10(pval_allChan),'ko');
            % Remove unqualified channels
            vNonQualifyingChan = find(~ismember(1:length(goodChan),vChan));
            plot(vNonQualifyingChan, log10(pval_allChan(vNonQualifyingChan)),'co');
            if TargFreq == 23
                plot(vChan(iMinPval_qualifying)*[1 1], ylim,'g--')
            else
                plot(vChan(iMinPval_qualifying)*[1 1], ylim,'r--')
            end
            %             if TargFreq == 23
            %             plot(iMin_pval*[1 1], ylim,'g--')
            %             else
            %             plot(iMin_pval*[1 1], ylim,'r--')
            %             end
            ylabel('log10(pval)','FontSize',12)
            xlabel('Channel','FontSize',12)
            legend('all pvals','nonQualityChan',['minPval:' num2str(vChan(iMinPval_qualifying))],'Location','SouthOutside','Orientation','Horizontal')
            title(['S1 ' Effect(iEffect,:) ' main effect,at ' num2str(TargFreq) 'Hz response'],'Fontsize',15);
            print(8, '-dpng', [Effect(iEffect,:) ' main effect response at ' num2str(TargFreq) 'Hz'])
        end
        
        %% After determining the matched experiment and corresponding channel number.
        % plot the power spectrum at this channel for this exp
        if fig9
            %% Define frequency
            load(fullfile(dir_sec, loadname(1).name));
            FreqPwr = data.freq{1};
            %% Plot power spectrum
            f9 = figure(9); clf
            set(f9,'Position',[1000 1000 800 1000])
            for iFile = 1:length(loadname)
                subplot(xlength,ylength,iFile)
                % plot(data.freq{1},mean(allData(iFile,:,:),3))
                plot(FreqPwr,mean(allData(iFile,:,:),3))
                %                 switch Effect(iEffect,:)
                %                     case {'  a'}
                FreqRange = [20 25];
                PwrRange = [0 35];
                
                %                     case {'  b'}
                %                         FreqRange = [195 205];
                %                         PwrRange = [5 15];
                %                     otherwise
                %                         % FreqRange = [0 250];
                %                         % PwrRange = [0 40];
                %                         FreqRange = [20 25];
                %                         PwrRange = [0 40];
                %                 end
                xlim(FreqRange)
                ylim(PwrRange)
                set(gca,'xtick',[])
                set(gca,'ytick',[])
            end
            set(gcf,'NextPlot','add','PaperPosition', [5 10 23 23]);
            axes;
            h = title([MatchFile ' Power Spectrum effect of ' Effect(iEffect,:)]);
            set(gca,'VIsible','off')
            set(h,'VIsible','on','interpret','none','VerticalAlignment','bottom')
            
            x=axes('units','normalized','position',[0.14 0.05 0.75 0],'xlim',[1 5],'color','none');
            xlabel(x,'200Hz stimulus strength')
            set(gca,'xtick',1:xlength,'XTickLabel',ConValue200)
            x2=axes('units','normalized','position',[0.14 0.095 0.75 0],'xtick',[],'color','b');
            xlabel(x2,['Freq Range is ' num2str(FreqRange(1)) ' to ' num2str(FreqRange(2))])
            
            y=axes('units','normalized','position',[0.07 0.09 0 0.8],'ylim',[1 5],'color','none');
            ylabel(y,'23Hz stimulus strength')
            set(gca,'ytick',1:ylength,'YTickLabel',flipud(ConValue23))
            y2=axes('units','normalized','position',[0.115 0.1 0 0.8],'ytick',[],'color','b');
            ylabel(y2,['Power Spectrum Strength from ' num2str(PwrRange(1)) ' to ' num2str(PwrRange(2))])
            print(f9,'-dpng',[MatchFile 'PwrSpec_effect_' Effect(iEffect,:)])
            
            %%
            %             %[tmp, iTmp] = min(abs(data.freq{1}-23))
            %             [tmp, iTmp] = min(abs(FreqPwr-TargFreq));
            %             f10 = figure(10); clf
            %             set(f10,'Position',[1000 1300 1000 1000])
            %             for iFile = 1:length(loadname)
            %                 subplot(xlength,ylength,iFile)
            %                 plot(squeeze(allData(iFile,iTmp,:)))
            %
            %                 PwrRange = [-10 40];
            %                 ylim(PwrRange)
            %                 set(gca,'xtick',[])
            %                 set(gca,'ytick',[])
            %             end
            %
            %             set(gcf,'NextPlot','add','PaperPosition', [10 10 20 20]);
            %             axes;
            %             h = title([MatchFile ' Power Spectrum effect of ' Effect(iEffect,:) ', response at ' num2str(FreqPwr(iTmp)) 'Hz']);
            %             set(gca,'VIsible','off')
            %             set(h,'VIsible','on','interpret','none','VerticalAlignment','bottom')
            %
            %             x=axes('units','normalized','position',[0.14 0.05 0.75 0],'xlim',[1 5],'color','none');
            %             xlabel(x,'200Hz stimulus strength')
            %             set(gca,'xtick',1:xlength,'XTickLabel',ConValue200)
            %             x2=axes('units','normalized','position',[0.14 0.095 0.75 0],'xtick',[],'color','b');
            %             xlabel(x2,['No. of Trials is ' num2str(size(allData(iFile,iTmp,:),3))])
            %
            %             y=axes('units','normalized','position',[0.07 0.09 0 0.8],'ylim',[1 5],'color','none');
            %             ylabel(y,'23Hz stimulus strength')
            %             set(gca,'ytick',1:ylength,'YTickLabel',flipud(ConValue23))
            %             y2=axes('units','normalized','position',[0.115 0.1 0 0.8],'ytick',[],'color','b');
            %             ylabel(y2,['Power Spectrum display from ' num2str(PwrRange(1)) ' to ' num2str(PwrRange(2))])
            %
            %             print(f10,'-dpng',[MatchFile 'PwrSpec_AllTrial_effect_' Effect(iEffect,:)])
            
            %% PLot 2D power spectrum across all conditions
            %             [tmp, iTmp] = min(abs(FreqPwr-TargFreq));
            %             f11 = figure(11); clf
            %             set(f11,'Position',[1000 1000 1000 1000])
            %             xticks = str2num(ConValue200);
            %             yticks = str2num(ConValue23);
            %             imagesc(1:xlength, 1:ylength,reshape(squeeze(mean(allData(:,iTmp,:),3)'),xlength,ylength)')
            %             xlabel('200 Hz stimulus strength')
            %             ylabel('23 Hz stimulus strength')
            %             set(gca,'xtick',1:xlength,'xticklabel',num2str(xticks),'ytick',1:ylength,'yticklabel',num2str(yticks))
            %             switch Effect(iEffect,:)
            %                 case {'  a'}
            %                     title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pa(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),2))]},'interpret','none','FontSize',13)
            %                 case {'  b'}
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pb(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),1))]},'interpret','none','FontSize',13)
            %                 case {'  c'}
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pc(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),3))]},'interpret','none','FontSize',13)
            %                 case {' ab'}
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pa(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),2))],...
            %                 ['minimum Pb(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),1))]},'interpret','none','FontSize',13)
            %                 case {' ac'}
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pa(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),2))],...
            %                 ['minimum Pc(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),3))]},'interpret','none','FontSize',13)
            %                 case {' bc'}
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pb(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),1))],...
            %                 ['minimum Pc(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),3))]},'interpret','none','FontSize',13)
            %                 otherwise
            %                      title({[MatchFile ' PwrSpec effect:' Effect(iEffect,:) ',response at ' num2str(FreqPwr(iTmp)) 'Hz'],...
            %                 ['minimum Pa(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),2))],...
            %                 ['minimum Pb(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),1))],...
            %                 ['minimum Pc(not log10): ' num2str(good_pval(iTargFreq,vChan(iMinPval_qualifying),3))]},'interpret','none','FontSize',13)
            %             end
            %             colorbar
            %
            %             print(f11,'-dpng',[MatchFile 'PwrSpec2D_AllCond_effect_' Effect(iEffect,:)])
            
            %             %% Plot 3D power spectrum strength in all conditions
            %             [tmp, iTmp] = min(abs(FreqPwr-TargFreq));
            %             f12 = figure(12); clf
            %             set(f12,'Position',[1000 1000 1000 1000])
            %             xticks = str2num(ConValue200);
            %             % xticks = [0 19 40 79 159]';
            %             yticks = str2num(ConValue23);
            %             % yticks = [0 2 4 7 15]';
            %             bar3(reshape(squeeze(mean(allData(:,iTmp,:),3)),length(xticks),length(yticks))')
            %             % bar3(reshape(squeeze(mean(allData(:,iTmp,:),3)),length(xticks),length(yticks))')
            %             xlabel('200 Hz stimulus strength')
            %             ylabel('23 Hz stimulus strength')
            %             zlabel('Power Spectrum')
            %             view(-135,30)
            %             set(gca,'xtick',1:xlength,'xticklabel',num2str(xticks),'ytick',1:ylength,'yticklabel',num2str(yticks))
            %             title([MatchFile ' PwrSpec response at ' num2str(FreqPwr(iTmp)) 'Hz'],'interpret','none')
            %
            %             print(f12,'-dpng',[MatchFile 'PwrSpec3D_AllCond_effect_' Effect(iEffect,:)])
            %
        end
        
    end
    %% Plot channel amount distribution
    if figPie
        figure(13);clf
        pie(ChanDistribution,{'23Hz','200Hz','Interaction','23Hz and 200Hz','23Hz and interaction','200Hz and interaction','All 3 factors'})
        title({[num2str(TargFreq) 'Hz response S2' ],'Distributions of channel amount satisfy main effect of:'},'interpret','none','FontSize',15)
        print(13,'-dpng','Distribution_channel_underEffect')
    end
end
%% Phoebe method to find the exp and real channel and plot

%         [a,b] = min(good_pval(iTargFreq,vChan,2));
%         iMinPvalChan = vAllChan(b);
%         [m,n] = find(all_pval == a);
%         FakeChan = mod(n,size(all_pval,2)); % after concatenating channel number is 1620
%         % Determine the channel number regardless of 3 dimention matrix
%         RealChan = mod(FakeChan,size(pvals,1));% Each data has 180 chan for S1 and 112 for S2
%         if RealChan == 0
%             RealChan = size(pvals,1);
%         end
%         RealExp = floor(FakeChan/size(pvals(:,:,2),1))+1;
%         disp(['The minimum pval belongs to experiment ' filename(RealExp,:) ',at channel:' num2str(RealChan)]);
%         if size(RealExp)~=1 % In case not only one datasets have the same minimum p value
%             disp('Multi RealExp Match')
%         end
%             FreqRange = [0 220];
%             MatchFile = [filename(RealExp,:) A{:}];
%             dir_sec = [data_dir(1:end-48) 'epoched_rsampsl_biprref_evkresp_cmtspwr/'];
%             loadname = dir(fullfile(dir_sec, [MatchFile '*.mat']));
%             % The first data (minimum both zero)
%             iFile = 1;
%             load(fullfile(dir_sec, loadname(iFile).name))
%             %     load(fullfile(dir_sec, loadname(iFile).name),'struct.trial','struct.freq')
%             [~, fstart_ind] = find_closest(data.freq{1}, FreqRange(1));
%             [~, fend_ind]   = find_closest(data.freq{1}, FreqRange(2));
%             Y(iFile,:,iEffect) = mean(data.trial(RealChan,fstart_ind:fend_ind, :),3);%define Y with 200Hz strength,23Hz no strength
%             X(iFile,:,iEffect) = mean(data.trial(RealChan,fstart_ind:fend_ind, :),3);%define X with 23Hz strength,200Hz no strength
%             % The rest of data
%             for iFile = 2:size(loadname,1)
%                 load(fullfile(dir_sec, loadname(iFile).name))
%                 [~, fstart_ind] = find_closest(data.freq{1}, FreqRange(1));
%                 [~, fend_ind]   = find_closest(data.freq{1}, FreqRange(2));
%                 switch loadname(1,1).name(29:31) %The strngth of 23Hz stim
%                     case {'000'}
%                         Y(iFile,:,iEffect) = mean(data.trial(RealChan,fstart_ind:fend_ind, :),3);
%                         %define Y with 200Hz strength,23Hz no strength
%                     otherwise
%                         X(iFile,:,iEffect) = mean(data.trial(RealChan,fstart_ind:fend_ind, :),3);
%                         %define X with 23Hz strength,200Hz no strength
%                 end
%             end
