%% init

close all; clc
clear; clear global

global id_ad id_noad id_isi2 id_isi3 id_ori % declare all global var for single dataset
global frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list
global root_path

root_path = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter';
% cd(root_path)

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

% caiman test data: grating gcamp6s 
% Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002
dataset_list = struct;
dataset_list.mouse = [1329]; 
dataset_list.date = [201209];
dataset_list.area = {'V1'};

%% load [xls, timecourse, stim]

% for iset = 1 %: length(dataset_list.date)
iset = 1
save_flag = 1; % toggle this to save/skip all .mat creation below

date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse];
area = dataset_list.area{1,iset}
[npSub_tc, frame_rate, input_behav, info] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);

areamousedate = [area '_' imouse '_' date];
result_folder = [root_path, areamousedate];
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

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
if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end
close

range_base = 1:3; range_resp = 9:11;
% prompt = 'base window = 1:3. what is resp window? '; range_resp = input(prompt); close

%% response to adapter & targets. get trace 
% dfof_ad = ncell x 1. dfof_tg = ncell x nori x nisi

[dfof_ad, dfof_ad_sem, dfof_ad_std] = dfof_resp(dfof_align_ad, 'ad', 0); % 0 to prevent saving dfof_ad vs dfof_tg separately
[dfof_tg, dfof_tg_sem, dfof_tg_std] = dfof_resp(dfof_align_tg, 'tg', 0);
if save_flag; save dfof.mat dfof_ad dfof_ad_sem dfof_ad_std dfof_tg dfof_tg_sem dfof_tg_std; end 

% trace = ncell x nori x nisi3 [noad 750 250]
[trace_avg, trace_sem] = trace_grand_avg(dfof_align_ad, save_flag);

%% visually driven cells
% cells responsive to ad / noad tg (all oris)

sig_alpha = 0.01;
[sig_vis_ad, p_vis_ad, ~] = vis_cell_criteria(dfof_align_ad, 'ad', sig_alpha);
[sig_vis_noad_tg, p_vis_noad_tg, ~] = vis_cell_criteria(dfof_align_tg, 'tg_any', sig_alpha);
vis_cell_ad = logical(sig_vis_ad');
vis_cell_noad_tg = logical(sum(sig_vis_noad_tg, 2));
vis_cell_any = logical(vis_cell_ad | vis_cell_noad_tg);
sum(vis_cell_ad) % vis driven by ad = 118 manual
sum(vis_cell_noad_tg) % vis driven by noad tg = 87 manual
sum(vis_cell_any) % vis driven by any = 123 manual

% find(vis_cell_ad==0) % not vis driven by ad
% find(vis_cell_noad_tg==0) % not vis driven by noad tg
% find(~vis_cell_ad & ~vis_cell_noad_tg) % not vis driven by anything

%% well-fit cells
% cells whose noad-tg 90% bootstraps are within 22.5 deg of all-trials-included fit

% bootstrap_file = fullfile(result_folder, 'fit_bootstrap.mat');
% if exist(bootstrap_file, 'file'); load(bootstrap_file, 'well_fit_cell')
% else
    cd(result_folder); nrun = 1000; save_flag = 1;
    well_fit_cell = well_fit_cell_criteria(dfof_align_tg, nrun, save_flag); 
% end
% sum(well_fit_cell)

%% fit tuning
% fit tuning under conditions = ncell x nparam x nisi [noad vs ad750 vs ad250]
[fit_param, ori_pref] = fit_tuning(dfof_tg, save_flag);

%% cell property

if save_flag
    save cell_property_loose.mat vis_cell_ad vis_cell_noad_tg sig_vis_ad sig_vis_noad_tg...
    ori_pref % well_fit_cell
end

% end