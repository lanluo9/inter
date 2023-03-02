%% init

close all; clc; clear
cd C:\Users\lan\Documents\repos\inter\mat

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

master_xls = readtable('C:\Users\lan\Documents\repos\inter\mat\adp_dataset_master.xlsx');

% dataset_list = struct;
% dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
% dataset_list.date = [200803, 200804, 200806,...
%                     200720, 200721, 200723, ...
%                     200728, 200729];
% dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};

%% load [xls, timecourse, stim]

% for iset = 1 %: length(dataset_list.date)
% iset
save_flag = 0; % toggle this to save/skip all .mat creation below

clear id_ad id_noad id_isi2 id_isi3 id_ori
clear frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list
global id_ori id_ad id_noad id_isi2 id_isi3
global frame_rate range_base range_resp ncell ntrial trial_len_min nori ori_list nisi

% date = num2str(dataset_list.date(iset))
% mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse];
% area = dataset_list.area{1,iset}
irow = find(master_xls.date == 210616); 
date = num2str(master_xls.date(irow));
assert(strcmp(date,'210616'))
mouse = master_xls.mouse(irow); imouse = ['i', num2str(mouse)];
area = master_xls.area(irow); area = area{1};

[npSub_tc, frame_rate, input_behav, info, result_folder] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);
cd(result_folder)

%% params & indexing trials
% index by adapter contrast, target ori, isi

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps 
% final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc);

adapter_id = cell2mat(input_behav.tstimOne); adapter_id = adapter_id(1:ntrial);
adapter_list = unique(adapter_id); n_adapter = length(adapter_list);
target_id = cell2mat(input_behav.tstimTwo); target_id = target_id(1:ntrial);
target_list = unique(target_id); n_target = length(target_list);
% verify randStim1_doSameStim2 protocol
assert(input_behav.doRandStimOne == 1 & input_behav.doSameStims == 1)
assert(sum(adapter_id == target_id) == length(target_id))

% contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast); 
% id_noad = find(contrast_ad == 0); id_ad = find(contrast_ad == 1); 
% id_noad(id_noad > ntrial) = []; id_ad(id_ad > ntrial) = []; 
id_noad = []; id_ad = 1:ntrial; 

frame_ad = double(cell2mat(input_behav.cStimOneOn)); frame_ad = frame_ad(1:ntrial);
frame_ad_off = double(cell2mat(input_behav.cStimOneOn)); frame_ad_off = frame_ad_off(1:ntrial);
frame_tg = double(celleqel2mat_padded(input_behav.cStimTwoOn)); frame_tg = frame_tg(1:ntrial);
isi_seq = frame_tg - frame_ad_off; 
trial_len_min = min(unique(diff(frame_ad)));

nisi = length(unique(frame_tg - frame_ad));
% id_750 = find(isi_seq > mean(isi_seq)); id_250 = find(isi_seq < mean(isi_seq)); 
% id_750(id_750 > ntrial) = []; id_250(id_250 > ntrial) = []; 
id_750 = []; id_250 = 1:ntrial; 
id_ad750 = intersect(id_ad, id_750); id_ad250 = intersect(id_ad, id_250);
id_isi2 = {id_ad750, id_ad250}; 
id_isi3 = {id_noad, id_ad750, id_ad250};

ori_seq = adapter_id;
ori_list = adapter_list; 
nori = length(ori_list); id_ori = cell(nori, 1);
for iori  = 1 : nori
    id_ori{iori} = find(ori_seq == ori_list(iori)); 
end

%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

tc_align_ad = align_tc(frame_ad, npSub_tc);
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad);
dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad);

%% set resp window
% find base window & resp window

t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); t_tg = squeeze(nanmean(t(:,:), 1)); 
range = 50; plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize')); legend('ad align', 'targ align')
% if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end
% close

range_base = 1:3; range_resp = 9:12;
% prompt = 'base window = 1:3. what is resp window? '; range_resp = input(prompt); close

%% response to adapter & targets. get trace 
% dfof_ad = ncell x nstim. dfof_tg = ncell x nstim

[dfof_ad, dfof_ad_sem, dfof_ad_std] = dfof_resp(dfof_align_ad, 'tg', 0); % tg mode aka separate diff stim images, but use adapter resp
[dfof_tg, dfof_tg_sem, dfof_tg_std] = dfof_resp(dfof_align_tg, 'tg', 0);
dfof_ad = dfof_ad(:,:,3); dfof_ad_sem = dfof_ad_sem(:,:,3); dfof_ad_std = dfof_ad_std(:,:,3); % keep isi=250 only
dfof_tg = dfof_tg(:,:,3); dfof_tg_sem = dfof_tg_sem(:,:,3); dfof_tg_std = dfof_tg_std(:,:,3);
if save_flag; save dfof.mat dfof_ad dfof_ad_sem dfof_ad_std dfof_tg dfof_tg_sem dfof_tg_std; end 

% trace = ncell x nori x nisi3 [noad 750 250]
[trace_avg, trace_sem] = trace_grand_avg(dfof_align_ad, 0);
trace_avg = squeeze(trace_avg(:,:,3,:)); trace_sem = squeeze(trace_sem(:,:,3,:));
if save_flag; save trace_aligned.mat trace_avg trace_sem; end

%% visually driven cells
% cells responsive to adapter aka stimOne, categorized by adapter identity

sig_alpha = 0.01;
% [sig_vis_ad, p_vis_ad, ~] = vis_cell_criteria(dfof_align_ad, 'ad', sig_alpha);
% vis_cell_ad = logical(sig_vis_ad');
id_noad = id_ad; % pretend every trial is noad bc vis_cell_criteria tg_any is for noad tg
[sig_vis_ad, p_vis_ad, ~] = vis_cell_criteria(dfof_align_ad, 'tg_any', sig_alpha);
vis_cell_ad = logical(sum(sig_vis_ad, 2));

sum(vis_cell_ad)/length(vis_cell_ad) % vis driven by any stim
sum(sig_vis_ad)/length(sig_vis_ad) % vis driven by each stim
histogram(sum(sig_vis_ad,2)) % how many neurons respond to 0...7 stimuli?

if save_flag
    save vis_driven.mat sig_vis_ad p_vis_ad vis_cell_ad
end

% find(vis_cell_ad==0) % not vis driven by ad
% find(vis_cell_noad_tg==0) % not vis driven by noad tg
% find(~vis_cell_ad & ~vis_cell_noad_tg) % not vis driven by anything

%% feature matrix 
% ntrial x ncell matrix
% side columns (img_id = 0-6, rep_num = 0 for ad | 1 for tg)

% feature_matrix
feature_ad = mean(dfof_align_ad(:,:,range_resp), 3) - mean(dfof_align_ad(:,:,range_base), 3);
feature_ad = feature_ad';
feature_tg = mean(dfof_align_tg(:,:,range_resp), 3) - mean(dfof_align_tg(:,:,range_base), 3);
feature_tg = feature_tg';

stim_id = double((adapter_id-1)');
rep_num_0 = zeros(size(stim_id)); % rep0 = adapter
rep_num_1 = ones(size(stim_id)); % rep1 = target

% feature matrix with side columns
dim_red_ad = [stim_id, rep_num_0, feature_ad];
dim_red_tg = [stim_id, rep_num_1, feature_tg];

adp_abs = mean((feature_tg - feature_ad) ./ feature_ad, 1); % absolute adaptation index
adp_abs = abs(adp_abs)'; % average adp over stimuli
histogram(adp_abs, 100)
adp_vis = [vis_cell_ad, adp_abs];

save dim_red.mat dim_red_ad feature_ad dim_red_tg feature_tg adp_vis

%%
% %% well-fit cells
% % cells whose noad-tg 90% bootstraps are within 22.5 deg of all-trials-included fit
% 
% % bootstrap_file = fullfile(result_folder, 'fit_bootstrap.mat');
% % if exist(bootstrap_file, 'file'); load(bootstrap_file, 'well_fit_cell')
% % else
%     cd(result_folder); nrun = 1000; save_flag = 1;
%     well_fit_cell = well_fit_cell_criteria(dfof_align_tg, nrun, save_flag); 
% % end
% % sum(well_fit_cell)

% %% fit tuning
% % fit tuning under conditions = ncell x nparam x nisi [noad vs ad750 vs ad250]
% [fit_param, ori_pref] = fit_tuning(dfof_tg, save_flag);

% %% cell property
% 
% if save_flag
%     save cell_property_loose.mat vis_cell_ad vis_cell_noad_tg sig_vis_ad sig_vis_noad_tg...
%     ori_pref well_fit_cell
% end

% end