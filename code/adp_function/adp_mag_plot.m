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

%% adaptation magnitude
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
    {'mean',@(adp_a0t0_merge)  prctile(adp_a0t0_merge,50)});

%% test necessity of thresholding (even after vis_cell_ad ensures dfof_ad)



%% tuning bias after adp | Jin2019 Fig 2F
% % y axis: change of distance of pref ori from adapter after adaptation (with_ad - no_ad)
% dis_pref = ori_pref(well_fit_cell,:); 
% dis_pref(dis_pref > 90) = abs(dis_pref(dis_pref > 90) - 180);
% dis_pref_change = dis_pref(:,2:3) - dis_pref(:,1);
% 
% % x axis: sort distance into 3 bins
% dis_pref_bin = dis_pref(:,1);
% dis_pref_bin(dis_pref_bin < 22.5) = 0; 
% dis_pref_bin(dis_pref_bin >= 22.5 & dis_pref_bin <= 67.5) = 45; 
% dis_pref_bin(dis_pref_bin > 67.5) = 90; 
% 
% % histogram(dis_pref_change((dis_pref_bin==0),2),53) % dist of 0 bin: not very pos (repulsive)
% 
% bin_list = unique(dis_pref_bin); nbin = length(unique(dis_pref_bin)); 
% x = []; y = [];
% for ibin = 1 : nbin
%     for iisi = 1 : nisi
%         x(ibin) = bin_list(ibin);
%         y(ibin, iisi) = mean(dis_pref_change(dis_pref_bin == x(ibin), iisi));
%         y2(ibin, iisi) = median(dis_pref_change(dis_pref_bin == x(ibin), iisi));
%     end
% end
