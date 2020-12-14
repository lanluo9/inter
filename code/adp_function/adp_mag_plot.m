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
cd C:\Users\lan\Documents\repos\inter\plot\

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
subplot(3,1,1)
adp_mean_750 = movmean(adp_sort_ad(:,2), rolling_win); plot(adp_mean_750, 'b'); hold on;
adp_mean_250 = movmean(adp_sort_ad(:,3), rolling_win); plot(adp_mean_250, 'r'); 
grid minor; yline(0, 'g'); xlabel('cells sorted by adapter response')
ylabel('adp moving mean'); ylim([-1,1])
subplot(3,1,2)
adp_std_750 = movstd(adp_sort_ad(:,2), rolling_win); plot(adp_std_750); hold on;
adp_std_250 = movstd(adp_sort_ad(:,3), rolling_win); plot(adp_std_250); grid minor; yline(1, 'g');
ylabel('adp moving std'); ylim([0,2])
subplot(3,1,3)
scatter(adp_sort_ad(:,1),adp_sort_ad(:,2),20,'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
scatter(adp_sort_ad(:,1),adp_sort_ad(:,3),20,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); grid minor;
xlabel('adapter response'); ylabel('adaptation index')
saveas(gcf, ['adp movmean movstd across dfof_ad & scatter'], 'jpg'); close

cutoff_index = 90; % or 50
dfof_ad_cutoff = adp_sort_ad(cutoff_index, 1);
disp(['movmean diverge bc isi & movstd decrease til <1. threshold determined: ' num2str(dfof_ad_cutoff)])

%% map scatter x axis onto mov mean std
rolling_win = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
adp_mean_750 = movmean(adp_sort_ad(:,2), rolling_win); plot(adp_sort_ad(:,1), adp_mean_750, 'b'); hold on;
adp_mean_250 = movmean(adp_sort_ad(:,3), rolling_win); plot(adp_sort_ad(:,1), adp_mean_250, 'r'); 
scatter(adp_sort_ad(:,1),adp_mean_750,10,'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
scatter(adp_sort_ad(:,1),adp_mean_250,10,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
grid minor; yline(0, 'g'); xlabel('adapter response')
ylabel('adp moving mean'); ylim([-1,1])
subplot(2,1,2)
adp_std_750 = movstd(adp_sort_ad(:,2), rolling_win); plot(adp_sort_ad(:,1), adp_std_750); hold on;
adp_std_250 = movstd(adp_sort_ad(:,3), rolling_win); plot(adp_sort_ad(:,1), adp_std_250); grid minor; yline(1, 'g');
ylabel('adp moving std'); ylim([0,2])
saveas(gcf, ['adp movmean movstd across dfof_ad all'], 'jpg'); close

zoomto_idx = find(adp_sort_ad(:,1) > 0.05, 1);
rolling_win = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
adp_mean_750 = movmean(adp_sort_ad(:,2), rolling_win); plot(adp_sort_ad(1:zoomto_idx,1), adp_mean_750(1:zoomto_idx), 'b'); hold on;
adp_mean_250 = movmean(adp_sort_ad(:,3), rolling_win); plot(adp_sort_ad(1:zoomto_idx,1), adp_mean_250(1:zoomto_idx), 'r'); 
scatter(adp_sort_ad(1:zoomto_idx,1),adp_mean_750(1:zoomto_idx),30,'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
scatter(adp_sort_ad(1:zoomto_idx,1),adp_mean_250(1:zoomto_idx),30,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
grid minor; yline(0, 'g'); xlabel('adapter response')
ylabel('adp moving mean'); ylim([-1,1])
subplot(2,1,2)
adp_std_750 = movstd(adp_sort_ad(:,2), rolling_win); plot(adp_sort_ad(1:zoomto_idx,1), adp_std_750(1:zoomto_idx)); hold on;
adp_std_250 = movstd(adp_sort_ad(:,3), rolling_win); plot(adp_sort_ad(1:zoomto_idx,1), adp_std_250(1:zoomto_idx)); grid minor; yline(1, 'g');
ylabel('adp moving std'); ylim([0,2])
saveas(gcf, ['adp movmean movstd across dfof_ad zoomin'], 'jpg'); close

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

dis_pref = []; area_merge = [];
for iset = 1 : nset
    well_fit_cell = set(iset).cell_property.well_fit_cell;
    ncell_set = sum(well_fit_cell); % qualified cell per set
    
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set,1) * areacode;
    area_merge = [area_merge; temp];
    
    temp = set(iset).cell_property.ori_pref(well_fit_cell, :);
    dis_pref = [dis_pref; temp];    
end

% y axis: change of distance of pref ori from adapter after adaptation (with_ad - no_ad)
dis_pref(dis_pref > 90) = abs(dis_pref(dis_pref > 90) - 180);
dis_pref_change = dis_pref(:,2:3) - dis_pref(:,1);
% x axis: sort distance into 3 bins
dis_pref_bin = dis_pref(:,1);
dis_pref_bin(dis_pref_bin < 22.5) = 0; dis_pref_bin(dis_pref_bin > 67.5) = 90; 
dis_pref_bin(dis_pref_bin >= 22.5 & dis_pref_bin <= 67.5) = 45; 

T_dis = table(area_merge, dis_pref_bin, dis_pref_change);
stat_dis = grpstats(T_dis,{'area_merge', 'dis_pref_bin'},{'mean','sem', 'min','max'},...
    'DataVars','dis_pref_change')

[stat_mean_dis,stat_median_dis]= grpstats(dis_pref_change,dis_pref_bin, ...
    {'mean',@(dis_pref_change)  prctile(dis_pref_change,50)})

%% plot for V1 only (low ncell for LM & LI)

dis_pref_change_v1 = dis_pref_change(area_merge==1,:);
dis_pref_bin_v1 = dis_pref_bin(area_merge==1,:);

T_dis_v1 = table(dis_pref_bin_v1, dis_pref_change_v1);
stat_dis = grpstats(T_dis_v1,{'dis_pref_bin_v1'},{'mean','sem', 'min','max'},...
    'DataVars','dis_pref_change_v1')

[stat_mean_dis_v1, stat_median_dis_v1]= grpstats(dis_pref_change_v1, dis_pref_bin_v1, ...
    {'mean',@(dis_pref_change_v1)  prctile(dis_pref_change_v1,50)})

%%

adp_area_avg = stat_dis.mean_dis_pref_change_v1; adp_area_sem = stat_dis.sem_dis_pref_change_v1;
ncell_area = stat_dis.GroupCount;
color_list = {[0,0,1], [1,0,0]}; clear h; figure
for iisi = 1:2
    h{iisi,1} = scatter(1:3, adp_area_avg(:, iisi), 5, color_list{iisi}, 'filled'); hold on
    errorbar(1:3, adp_area_avg(:, iisi), adp_area_sem(:, iisi), 'color', color_list{iisi}, 'LineStyle','none')
end

ylim([-10,10]); 
yl = ylim;
for itext = 1:3
    text(itext, yl(1) + 1, ...
        ['n=', num2str(ncell_area(itext))], 'HorizontalAlignment', 'center')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xlim([0.5, 3.5]); xticks(1:3); xticklabels({'0', '45', '90'}); 
xlabel('|pref - adapter|')
ylabel('delta pref ori');
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250', 'Location','northeast'); legend boxoff
% saveas(gcf, ['distance of pref ori from adapter bef or after adaptation'], 'jpg'); close

%% LM & LI less ori-tuned
% proportion of well_fit_cell in vis_cell_noad_tg

well_fit_merge = []; well_fit_check = []; area_merge = [];
for iset = 1 : nset
    vis_noad_tg(iset) = sum(set(iset).cell_property.vis_cell_noad_tg);
    areacode(iset) = dataset_list.areacode{iset};
    
    temp = set(iset).cell_property.well_fit_cell(set(iset).cell_property.vis_cell_noad_tg);
    well_fit(iset) = sum(temp);   
end
well_fit_ratio = well_fit ./ vis_noad_tg;

T = table(areacode', well_fit_ratio');
stat_adp = grpstats(T,{'Var1'},{'mean','sem','max','min'},'DataVars','Var2')

[stat_mean,stat_median]= grpstats(well_fit_ratio,areacode, {'mean',@(well_fit_ratio)  prctile(well_fit_ratio,50)})
% interestingly, proportion of vis_ad cells in all cells
% & proportion of well-fit cells in vis_noad_tg cells 
% both double down as we proceed from V1 to LM to LI

%% OSI? SNR?
