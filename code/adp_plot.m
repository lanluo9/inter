%% set up 

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\code

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806,...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
for iset = 1 : length(dataset_list.date)
    areamouse_seq{iset} = [dataset_list.area{1,iset} '_' num2str(dataset_list.mouse(iset))];
end

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
set = cell(nset, 4);
for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset))
    imouse = ['i', mouse];
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];

    result_prefix = 'C:\Users\lan\Documents\repos\inter\code\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
    cd(result_folder{iset});

    mat_files = dir('*.mat'); 
    for imat = 1:length(mat_files) 
        set{iset, imat} = load(mat_files(imat).name); 
    end
end

%% adaptation index

ndelta = 8;
ngap = 2;

adp = struct;
for iset = 1 : nset
    
% with-adapter / no-adapter resp to same targ with same isi
cell_list_now = find(set{iset, 1}.vis_driven); 
adp_ratio = zeros(length(cell_list_now), ndelta, ngap);

for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    for idelta = 1:ndelta
    for igap = 1:ngap
        dfof_equiv_ad = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_noad = set{iset, 2}.dfof_avg_merge(icell, idelta, 1);
        adp_ratio(ii, idelta, igap) = dfof_equiv_ad / dfof_equiv_noad - 1;
    end
    end
end

% with-ad targ0 vs no-ad targ0
cell_list_now = set{iset, 1}.pref_0_cell;
adp_ratio_targ0 = zeros(length(cell_list_now), ngap);

for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8; % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_noad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, 1);
%         dfof_equiv_noad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 9:11)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 1:3)),2);
        adp_ratio_targ0(ii, igap) = dfof_equiv_ad_targ / dfof_equiv_noad_targ - 1;
    end
end

% with-ad targ0 vs its own ad0
adp_ratio_a0t0 = zeros(length(cell_list_now), ngap);
load(fullfile(result_folder{iset}, 'pre-processing', 'resp_ad.mat'))
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8 % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_ad = dfof_avg_ad(icell, idelta);
        adp_ratio_a0t0(ii, igap) = dfof_equiv_ad_targ / dfof_equiv_ad - 1;
    end
end

adp(iset).adp_ratio = adp_ratio;
adp(iset).adp_ratio_targ0 = adp_ratio_targ0;
adp(iset).adp_ratio_a0t0 = adp_ratio_a0t0;

end



