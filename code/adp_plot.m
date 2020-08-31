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
 
set = cell(length(dataset_list.date), 4);
for iset = 1 : length(dataset_list.date)
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset))
    imouse = ['i', mouse];
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];

    result_prefix = 'C:\Users\lan\Documents\repos\inter\code\';
    result_folder = fullfile(result_prefix, areamousedate);
    cd(result_folder);

    mat_files = dir('*.mat'); 
    for imat = 1:length(mat_files) 
        set{iset, imat} = load(mat_files(imat).name); 
    end
end

%% adaptation index
% with-adapter / no-adapter resp to same targ with same isi

ndelta = 8;
ngap = 2;

iset = 1;
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

%%
all = sum(adp_ratio(:)>-Inf)
outlier_facil = find(adp_ratio > mean(adp_ratio(:) + 3*std(adp_ratio(:))) | adp_ratio < mean(adp_ratio(:) - 3*std(adp_ratio(:))))
facilitated = sum(adp_ratio(:)>0)
inhibited = sum(adp_ratio(:)<=0)

adp_ratio(outlier_facil) = NaN;
nanmean(adp_ratio(:))
nanstd(adp_ratio(:))

% t750 = squeeze(adp_ratio(:,:,1));
% t250 = squeeze(adp_ratio(:,:,2));
% errorbar(nanmean(t750, 1), nanstd(t750, 1), 'b'); hold on;
% errorbar(nanmean(t250, 1), nanstd(t250, 1), 'r');
% line([0,9], [0, 0], 'Color', 'g', 'LineWidth', 1);
% xlabel('ori');ylabel('adp ratio')

figure
% subplot(1,2,1)
histogram(adp_ratio(:), 30)
yl = ylim;
line([0, 0], [0, yl(2)], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
xlabel('adp index')
ylabel('count across cell ori isi')

%%
% check only with-ad targ0 vs ad || with-ad targ0 vs no-ad targ0
cell_list_now = pref_0_cell;

adp_ratio_targ0 = zeros(length(cell_list_now), ngap);
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
% for icell = 1:ncell
    for idelta = 8 % targ0 only! adp is ori-specific
        id_delta = find(delta_seq == delta_list(idelta));
        for igap = 1:ngap
            idx_now_ad_targ = intersect(intersect(id_gaps{igap}, id_delta), id_ad);
%             idx_now_noad_targ = intersect(id_delta, id_noad);
            dfof_equiv_ad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_ad_targ, range_resp)),2)...
                               - mean(squeeze(tc_trial_align_targ(icell, idx_now_ad_targ, range_base)),2);
            dfof_equiv_ad = mean(squeeze(tc_trial_align_ad(icell, idx_now_ad_targ, range_resp)),2)...
                          - mean(squeeze(tc_trial_align_ad(icell, idx_now_ad_targ, range_base)),2);
%             dfof_equiv_noad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, range_resp)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, range_base)),2);
%             adp_ratio_targ0(ii, igap) = mean(dfof_equiv_ad_targ(dfof_equiv_ad_targ>0)) / mean(dfof_equiv_ad(dfof_equiv_ad>0)) - 1;
            adp_ratio_targ0(ii, igap) = mean(dfof_equiv_ad_targ) / mean(dfof_equiv_ad) - 1;
        end
    end
end

all = sum(adp_ratio_targ0(:)>-Inf)
outlier_facil0 = find(adp_ratio_targ0 > mean(adp_ratio_targ0(:) + 3*std(adp_ratio_targ0(:))))
facilitated = sum(adp_ratio_targ0(:)>0)
inhibited = sum(adp_ratio_targ0(:)<=0)

adp_ratio_targ0(outlier_facil0) = NaN;
nanmean(adp_ratio_targ0, 1)
nanstd(adp_ratio_targ0, 1)

subplot(1,2,2)
histogram(adp_ratio_targ0(:), 10)
yl = ylim;
line([0, 0], [0, yl(2)], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
xlabel('adp index')
ylabel('count across cell ori isi')

