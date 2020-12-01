%% set up dir

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\mat

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806, ...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
dataset_list.areacode = {1,2,3, 1,2,3, 1,2};

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse]
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];
    result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
end

%% load 

set = struct;
for iset = 1 : nset
    cd(result_folder{iset});
    set(iset).dfof = load('dfof.mat');
    set(iset).cell_property = load('cell_property.mat');
    set(iset).trace_aligned = load('trace_aligned.mat');    
end

%% adaptation magnitude (using vis_cell_ad before thresholding)
% adp_mag = dfof_withad_tg0 vs dfof_ad

dfof_wad_tg0_merge = []; dfof_ad_merge = []; area_merge = [];
for iset = 1 : nset
    vis_cell_ad = set(iset).cell_property.vis_cell_ad;
    ncell_set = sum(vis_cell_ad); % qualified cell per set
    
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set,1) * areacode;
    area_merge = [area_merge; temp];
    
    temp = set(iset).dfof.dfof_ad(vis_cell_ad);
    dfof_ad_merge = [dfof_ad_merge; temp];    
    temp = squeeze(set(iset).dfof.dfof_tg(vis_cell_ad, 1, 2:3)); % target ori=0, with ad
    dfof_wad_tg0_merge = [dfof_wad_tg0_merge; temp];
end
adp_a0t0_merge = dfof_wad_tg0_merge ./ dfof_ad_merge - 1;

T = table(area_merge, adp_a0t0_merge);
stat_adp = grpstats(T,{'area_merge'},{'mean','sem', 'min','max'},'DataVars','adp_a0t0_merge')

[stat_mean,stat_median]= grpstats(adp_a0t0_merge,area_merge, ...
    {'mean',@(adp_a0t0_merge)  prctile(adp_a0t0_merge,50)})

%% test necessity of thresholding (even after vis_cell_ad ensures dfof_ad)
% moving adp_std as dfof_ad changes

adp_sort_ad = [dfof_ad_merge, adp_a0t0_merge];
adp_sort_ad = sortrows(adp_sort_ad, 1);
% scatter(adp_sort_ad(:,1), adp_sort_ad(:,2)); hold on;
% scatter(adp_sort_ad(:,1), adp_sort_ad(:,3));

rolling_win = 20;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
adp_mean_750 = movmean(adp_sort_ad(:,2), rolling_win); plot(adp_mean_750); hold on;
adp_mean_250 = movmean(adp_sort_ad(:,3), rolling_win); plot(adp_mean_250); grid minor; yline(0, 'g');
title('adp moving mean'); ylim([-1,1])
subplot(2,1,2)
adp_std_750 = movstd(adp_sort_ad(:,2), rolling_win); plot(adp_std_750); hold on;
adp_std_250 = movstd(adp_sort_ad(:,3), rolling_win); plot(adp_std_250); grid minor; yline(1, 'g');
title('adp moving std'); ylim([0,2])
% saveas(gcf, ['adp movmean movstd across dfof_ad'], 'jpg'); close

cutoff_index = 90; % or 50
dfof_ad_cutoff = adp_sort_ad(cutoff_index, 1);
disp(['movmean diverge bc isi & movstd decrease til <1. threshold determined: ' num2str(dfof_ad_cutoff)])

%% adp after thresholding

area_co = area_merge(dfof_ad_merge >= dfof_ad_cutoff);
adp_a0t0_co = adp_a0t0_merge((dfof_ad_merge >= dfof_ad_cutoff), :);
T_co = table(area_co, adp_a0t0_co);
stat_adp_co = grpstats(T_co,{'area_co'},{'mean','sem', 'min','max'},'DataVars','adp_a0t0_co')

[stat_mean_co,stat_median_co]= grpstats(adp_a0t0_co, area_co, ...
    {'mean',@(adp_a0t0_co)  prctile(adp_a0t0_co,50)})

%% plot

adp_area_avg = stat_adp_co.mean_adp_a0t0_co; adp_area_sem = stat_adp_co.sem_adp_a0t0_co;
ncell_area_co = stat_adp_co.GroupCount;
color_list = {[0,0,1], [1,0,0]}; clear h; figure
for iisi = 1:2
    h{iisi,1} = scatter(1:3, adp_area_avg(:, iisi), 5, color_list{iisi}, 'filled'); hold on
    errorbar(1:3, adp_area_avg(:, iisi), adp_area_sem(:, iisi), 'color', color_list{iisi}, 'LineStyle','none')
end

adp_area_avg = stat_adp.mean_adp_a0t0_merge; adp_area_sem = stat_adp.sem_adp_a0t0_merge;
ncell_area = stat_adp.GroupCount;
color_list = {[0,0.5,0.5], [0.5,0.5,0]}; 
for iisi = 1:2
    scatter(1:3, adp_area_avg(:, iisi), 5, color_list{iisi}, 'filled'); hold on
    errorbar(1:3, adp_area_avg(:, iisi), adp_area_sem(:, iisi), 'color', color_list{iisi}, 'LineStyle','none')
end

ylim([-1.2,1.2]); yl = ylim;
for itext = 1:3
    text(itext, yl(1) + 0.1, ...
        ['n=',num2str(ncell_area_co(itext)), '|', num2str(ncell_area(itext))], ...
        'HorizontalAlignment', 'center')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xlim([0.5, 3.5]); xticks(1:3); xticklabels({'V1', 'LM', 'LI'}); ylabel('adaptation index');
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250', 'Location','northeast'); legend boxoff
% saveas(gcf, ['adp mag bef or after more thresholding dfof_ad'], 'jpg'); close

%% tuning bias after adp | Jin2019 Fig 2F

% dfof_wad_tg0_merge = []; dfof_ad_merge = []; area_merge = [];
for iset = 1 : nset
    well_fit_cell = set(iset).cell_property.well_fit_cell;
    ncell_set = sum(well_fit_cell); % qualified cell per set
    
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set,1) * areacode;
    area_merge = [area_merge; temp];
    
    temp = set(iset).dfof.dfof_ad(vis_cell_ad);
    dfof_ad_merge = [dfof_ad_merge; temp];    
    temp = squeeze(set(iset).dfof.dfof_tg(vis_cell_ad, 1, 2:3)); % target ori=0, with ad
    dfof_wad_tg0_merge = [dfof_wad_tg0_merge; temp];
end

% y axis: change of distance of pref ori from adapter after adaptation (with_ad - no_ad)
dis_pref = ori_pref(well_fit_cell,:); 
dis_pref(dis_pref > 90) = abs(dis_pref(dis_pref > 90) - 180);
dis_pref_change = dis_pref(:,2:3) - dis_pref(:,1);

% x axis: sort distance into 3 bins
dis_pref_bin = dis_pref(:,1);
dis_pref_bin(dis_pref_bin < 22.5) = 0; 
dis_pref_bin(dis_pref_bin >= 22.5 & dis_pref_bin <= 67.5) = 45; 
dis_pref_bin(dis_pref_bin > 67.5) = 90; 

% histogram(dis_pref_change((dis_pref_bin==0),2),53) % dist of 0 bin: not very pos (repulsive)

bin_list = unique(dis_pref_bin); nbin = length(unique(dis_pref_bin)); 
x = []; y = [];
for ibin = 1 : nbin
    for iisi = 1 : nisi
        x(ibin) = bin_list(ibin);
        y(ibin, iisi) = mean(dis_pref_change(dis_pref_bin == x(ibin), iisi));
        y2(ibin, iisi) = median(dis_pref_change(dis_pref_bin == x(ibin), iisi));
    end
end
