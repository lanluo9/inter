%% init
close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\mat

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806,...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};

%% load [xls, timecourse, stim]

% for iset = 1 : length(dataset_list.date)
iset = 1
date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse];
area = dataset_list.area{1,iset}

[npSub_tc, frame_rate, input_behav, info, result_folder] = ...
    load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);

%% params & indexing trials

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps % final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc);

contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast); unique(contrast_ad); 
ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg);
ori_seq(ori_seq == 180) = 0;
ori_list = unique(ori_seq); nori = length(ori_list); 

frame_ad = double(cell2mat(input_behav.cStimOn)); frame_ad_off = double(cell2mat(input_behav.cStimOff));
frame_tg = celleqel2mat_padded(input_behav.cTargetOn); frame_tg = double(frame_tg);
isi_seq = frame_tg - frame_ad_off; 
nisi = length(unique(frame_tg - frame_ad));

id_750 = find(isi_seq > mean(isi_seq)); id_750(id_750 > ntrial) = []; id_250 = find(isi_seq < mean(isi_seq)); id_250(id_250 > ntrial) = [];
id_gaps = {id_750, id_250}; 
id_noad = find(contrast_ad == 0); id_noad(id_noad > ntrial) = []; id_ad = find(contrast_ad == 1); id_ad(id_ad > ntrial) = [];

%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

tc_align_ad = align_tc(frame_ad, npSub_tc, ncell, ntrial);
tc_align_tg = align_tc(frame_tg, npSub_tc, ncell, ntrial);
dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad, frame_rate, ncell, ntrial);
dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad, frame_rate, ncell, ntrial);

%% set resp window
% find base window & resp window

range = 50;
t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); t_tg = squeeze(nanmean(t(:,:), 1)); 
plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize')); legend('ad align', 'targ align')
% saveas(gcf, 'dfof align zoomin', 'jpg'); close

range_base = [1:3]; range_resp = [9:12];
% prompt = 'base window = 1:3. what is resp window? ';
% range_resp = input(prompt); close

