%% set up 

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

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
set = struct;

for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse]
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];
    result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
    cd(result_folder{iset});

    set(iset).dfof_ad = load('dfof_ad.mat', 'dfof_ad');
    set(iset).dfof_tg = load('dfof_tg.mat', 'dfof_tg');
    set(iset).vis_cell_ad = load('cell_property.mat', 'vis_cell_ad');
end

%% adaptation magnitude
% adp_mag = dfof_withad_tg0 vs dfof_ad
% or dfof_withad_tg0 vs dfof_noad_tg0


