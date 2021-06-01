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
    set(iset).trace = load('trace_aligned.mat');
    
    set(iset).cell_property = load('cell_property_loose.mat');
    set(iset).fit_bootstrap = load('fit_bootstrap_90perc.mat');
    set(iset).fit_tuning = load('fit_tuning_isi3.mat');
end
cd C:\Users\lan\Documents\repos\inter\plot\


%% 

area_merge = []; vis = []; well_fit = []; ori_perc_all = [];
R2 = []; SSE = []; sharp = [];
ori_pref = []; ori_orth = []; dfof_pref = []; dfof_orth = []; dfof_pref_avg = []; dfof_pref_std = [];
ori_round_to = 0 : 22.5 : 157.5;

for iset = 1 : nset
    well_fit_cell = set(iset).cell_property.well_fit_cell;
    well_fit = [well_fit; well_fit_cell];
    
    vis_tg_cell = sum(set(iset).cell_property.sig_vis_noad_tg,2);
    vis_ad_cell = set(iset).cell_property.sig_vis_ad';
    vis_cell = vis_tg_cell | vis_ad_cell;
    vis = [vis; vis_cell];
    
    ncell_set(iset) = length(well_fit_cell); 
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set(iset),1) * areacode;
    area_merge = [area_merge; temp];
    
    tt = set(iset).fit_bootstrap.ori_perc(:, 1);
%     length(tt)
    ori_perc_all = [ori_perc_all; tt];
    tt = set(iset).fit_tuning.fit_param(:, end, 1);
    R2 = [R2; tt];
    tt = set(iset).fit_tuning.fit_param(:, 3, 1);
    sharp = [sharp; tt];
    tt = set(iset).fit_tuning.fit_param(:, end-1, 1);
    SSE = [SSE; tt];
    
    temp = set(iset).cell_property.ori_pref(:, 1);
    temp = interp1(ori_round_to,ori_round_to,temp, 'nearest','extrap');
    ori_pref = [ori_pref; temp];
    temp2 = temp + 90; temp2(temp2>=180) = temp2(temp2>=180) - 180;
    ori_orth = [ori_orth; temp2];
    
    for icell = 1 : ncell_set(iset)
%         cell_idx = max(find(well_fit_cell, icell));
        ori_pref_idx = find(ori_round_to == temp(icell));
        ori_orth_idx = find(ori_round_to == temp2(icell));
        dfof_pref_cell = set(iset).dfof.dfof_tg(icell, ori_pref_idx, 1);
        dfof_orth_cell = set(iset).dfof.dfof_tg(icell, ori_orth_idx, 1);
        dfof_pref = [dfof_pref; dfof_pref_cell];
        dfof_orth = [dfof_orth; dfof_orth_cell];
    end
    
    for icell = 1 : ncell_set(iset)
%         cell_idx = max(find(well_fit_cell, icell));
        ori_pref_idx = find(ori_round_to == temp(icell));
        dfof_pref_cell = set(iset).dfof.dfof_tg(icell, ori_pref_idx, 1);
        dfof_pref_avg = [dfof_pref_avg; dfof_pref_cell];
        dfof_pref_std_cell = set(iset).dfof.dfof_tg_std(icell, ori_pref_idx, 1);
        dfof_pref_std = [dfof_pref_std; dfof_pref_std_cell];
    end
end
OSI = (dfof_pref - dfof_orth) ./ (dfof_pref + dfof_orth); % OSI = diff(pref, ortho) / sum(pref, ortho) 
coeff_var = dfof_pref_std ./ dfof_pref_avg;

%%
ori_perc_mat = 'C:\Users\lan\Documents\repos\inter\plot\CV SNR OSI R2 ori_perc by area - why HVA lack well fit\ori_perc_all_area.mat'
load(ori_perc_mat)

save corr_well_fit_HVA.mat area_merge vis well_fit ori_perc_all R2 SSE sharp OSI coeff_var
