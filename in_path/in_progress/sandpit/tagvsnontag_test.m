%% Load from file
load('/media/phoebeyou/My Passport/Spencers_Cat_Data/functions__/aux_files/compiled_anova_data')
    


%% find the indices of the non-tagged frequencies
[freq_keep, keep_i, hars] = no_resp_indices([23, 200], [], data(1).freq, 'logical');
% data(1).freq == data(2).freq





%% seperate data

for cArea = 1:2
    data(cArea).nontagged.pvals = data(cArea).good.pvals(:, keep_i, :);
    data(cArea).nontagged.freq = freq_keep;
    
    data(cArea).tagged.pvals = data(cArea).good.pvals(:, ~keep_i, :);
    data(cArea).tagged.freq = data(1).freq(~keep_i);
    
end


 save('compiled_anova_data_sepbytag', 'data')




%% analysis
q = 0.05;

% pID or pN?????????

for cArea = 1:2
    % FDR correction
    [pID,pN] = FDR(data(cArea).nontagged.pvals,q);
    
    p_lim = pID;
    
    % disp(['Non-tagged (S' num2str(cArea) '): ' num2str(pID)])
    
    
    % classify the channels
    p_thresh = zeros(size(data(cArea).nontagged.pvals));
    p_thresh(data(cArea).nontagged.pvals<p_lim) = 1;
    % convert binary to decimal - int/(:, :, 3) is LSB
    p_class = p_thresh;
    p_class(:, :, 2) = p_class(:, :, 2)*2;
    p_class(:, :, 1) = p_class(:, :, 1)*4;
    p_class = sum(p_class, 3);
    
    
    nontagged_p_class{cArea} = p_class;
    nontagged_p_lim(cArea) = p_lim;
    
    
    [pID,pN] = FDR(data(cArea).tagged.pvals,q);
    
    % disp(['Tagged (S' num2str(cArea) '): ' num2str(pID)])
    
    
    p_lim = pID;
    
    
    % classify the channels
    p_thresh = zeros(size(data(cArea).tagged.pvals));
    p_thresh(data(cArea).tagged.pvals<p_lim) = 1;
    % convert binary to decimal - int/(:, :, 3) is LSB
    p_class = p_thresh;
    p_class(:, :, 2) = p_class(:, :, 2)*2;
    p_class(:, :, 1) = p_class(:, :, 1)*4;
    p_class = sum(p_class, 3);
    
    tagged_p_class{cArea} = p_class;
    tagged_p_lim(cArea) = p_lim;
end

data(1).nontagged.thresh = nontagged_p_class{1};
data(2).nontagged.thresh = nontagged_p_class{2};
data(1).tagged.thresh    = tagged_p_class{1};
data(2).tagged.thresh    = tagged_p_class{2};

data(1).nontagged.plim = nontagged_p_lim(1);
data(2).nontagged.plim = nontagged_p_lim(2);
data(1).tagged.plim    = tagged_p_lim(1);
data(2).tagged.plim    = tagged_p_lim(2);


save('compiled_anova_data_sepbytag_threshold', 'data')

%%
figure(1)

imagesc(tagged_p_class{1})
fnum = [4, 11, 18, 24, 29, 36, 43, 50, 57, 64, 70, 76, 82, 89, 96, 103, 110, 116, 121, 128, 135, 142];


set(gca, 'XTick', fnum)
set(gca, 'XTickLabel', hars)

title({'Classification of tagged in S1', ['p < ' num2str(tagged_p_lim(1))]})

colorbar
tab_base = de2bi(0:7, 'left-msb');
colormap(tab_base)

frow = '23';
fcol = '200';
colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
        [fcol '+int'], [frow '+' fcol], 'all'})


print(gcf, '-dpng', 'tagged_S1_classifications')

%%
figure(2)

imagesc(tagged_p_class{2})
fnum = [4, 11, 18, 24, 29, 36, 43, 50, 57, 64, 70, 76, 82, 89, 96, 103, 110, 116, 121, 128, 135, 142];


set(gca, 'XTick', fnum)
set(gca, 'XTickLabel', hars)

title({'Classification of tagged in S2', ['p < ' num2str(tagged_p_lim(2))]})

colorbar
tab_base = de2bi(0:7, 'left-msb');
colormap(tab_base)

frow = '23';
fcol = '200';
colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
        [fcol '+int'], [frow '+' fcol], 'all'})

print(gcf, '-dpng', 'tagged_S2_classifications')

%%
figure(3)

imagesc(data(1).nontagged.freq, [], nontagged_p_class{1})


title({'Classification of nontagged in S1', ['p < ' num2str(nontagged_p_lim(1))]})

colorbar
tab_base = de2bi(0:7, 'left-msb');
colormap(tab_base)

frow = '23';
fcol = '200';
colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
        [fcol '+int'], [frow '+' fcol], 'all'})


print(gcf, '-dpng', 'nontagged_S1_classifications')

%%
figure(4)

imagesc(data(2).nontagged.freq, [], nontagged_p_class{2})



title({'Classification of nontagged in S2', ['p < ' num2str(nontagged_p_lim(2))]})

colorbar
tab_base = de2bi(0:7, 'left-msb');
colormap(tab_base)

frow = '23';
fcol = '200';
colorbar('YTickLabel', {'none', 'interaction', frow, [frow '+int'], fcol, ...
        [fcol '+int'], [frow '+' fcol], 'all'})


print(gcf, '-dpng', 'nontagged_S2_classifications')


