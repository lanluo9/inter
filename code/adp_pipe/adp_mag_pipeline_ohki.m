%% init

close all; clc; clear
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

for iset = 1 : length(dataset_list.date)
iset
save_flag = 0; % toggle this to save/skip all .mat creation below

clear id_ad id_noad id_isi2 id_isi3 id_ori
clear frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list
global id_ad id_noad id_isi2 id_isi3 id_ori % declare all global var for single dataset
global frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list

date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse];
area = dataset_list.area{1,iset}
[npSub_tc, frame_rate, input_behav, info, result_folder] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);
result_folder = [result_folder, '_ohki'];
if ~exist(result_folder, 'dir'); mkdir(result_folder); end; cd(result_folder)

%% params & indexing trials
% index by adapter contrast, target ori, isi

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps 
% final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc);

contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast); 
id_noad = find(contrast_ad == 0); id_ad = find(contrast_ad == 1); 
id_noad(id_noad > ntrial) = []; id_ad(id_ad > ntrial) = []; 

frame_ad = double(cell2mat(input_behav.cStimOn)); frame_ad_off = double(cell2mat(input_behav.cStimOff));
frame_tg = celleqel2mat_padded(input_behav.cTargetOn); frame_tg = double(frame_tg);
isi_seq = frame_tg - frame_ad_off; 
nisi = length(unique(frame_tg - frame_ad));
id_750 = find(isi_seq > mean(isi_seq)); id_250 = find(isi_seq < mean(isi_seq)); 
id_750(id_750 > ntrial) = []; id_250(id_250 > ntrial) = []; 
id_ad750 = intersect(id_ad, id_750); id_ad250 = intersect(id_ad, id_250);
id_isi2 = {id_ad750, id_ad250}; 
id_isi3 = {id_noad, id_ad750, id_ad250};

trial_len_min = min(unique(diff(frame_ad)));

ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg); ori_seq(ori_seq == 180) = 0;
ori_seq(end) = [];
ori_list = unique(ori_seq); 
nori = length(ori_list); id_ori = cell(nori, 1);
for iori  = 1 : nori
    id_ori{iori} = find(ori_seq == ori_list(iori)); 
end

%%

ori_seq_noad = ori_seq(id_noad);
ori_index = ori_seq_noad / 22.5 + 1;
ntarget_noad = length(ori_seq_noad);

pathname = 'C:\Users\lan\Documents\repos\inter\code\borrowed\GaborWavelet\res\';
load([pathname, 'F_stim.mat'])
nfeature = size(F_stim,1);
F_trial = zeros(nfeature, ntarget_noad);
for i = 1 : ntarget_noad
    F_trial(:,i) = F_stim(:,ori_index(i));
end
save feature_trial.mat F_trial

%%
%{
%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

tc_align_ad = align_tc(frame_ad, npSub_tc);
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = dfof_by_trial_base_ohki(tc_align_ad, npSub_tc, frame_ad);
dfof_align_tg = dfof_by_trial_base_ohki(tc_align_tg, npSub_tc, frame_ad);

%% set resp window
% find base window & resp window

% t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); t_ad = squeeze(nanmean(t(:,:), 1)); 
% t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); t_tg = squeeze(nanmean(t(:,:), 1)); 
% range = 50; plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
% grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize')); legend('ad align', 'targ align')
% if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end; close

range_base = 1:3; range_resp = 9:12;
% prompt = 'base window = 1:3. what is resp window? '; range_resp = input(prompt); close

%% response to 8 ori targets. visually-driven & ori-driven cell
% dfof_ad = ncell x 1. dfof_tg = ncell x nori x nisi

sig_alpha = 0.01; resp_threshold = 0.1;
[dfof_tg, dfof_tg_sem, dfof_tg_std, ori_driven, vis_driven] = ...
    dfof_tg_vis_ori_driven_ohki(dfof_align_tg, sig_alpha, resp_threshold);

% subplot(1,2,1); imagesc(ori_driven)
% subplot(1,2,2); imagesc(ori_driven(vis_driven>0, :))
% t = sum(ori_driven,1); bar(t)

vis_driven_cell = vis_driven>0;
sum(vis_driven_cell)/length(vis_driven_cell)
ori_driven_cell = vis_driven>0 & sum(ori_driven,2)>0;
sum(ori_driven_cell)/length(ori_driven_cell)

if save_flag; save dfof_tg.mat dfof_tg dfof_tg_sem dfof_tg_std; ...
    save cell_property.mat ori_driven ori_driven_cell vis_driven vis_driven_cell; end 

%% response to no-ad targets, z-scored for each cell

save_flag = 1;
[dfof_tg_noad, z_score] = dfof_tg_noad_zscore_ohki(dfof_align_tg, save_flag);

%}

end