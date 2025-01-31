%% init

close all; clc; 
clear; clear global

%% 

database_path = 'C:\Users\ll357\Documents\inter\';
master_xls = [database_path, 'data/mix50_grat1.csv'];
dataset_meta = readtable(master_xls);

% segmented_id = (dataset_meta.cellpose_seg == 1) | (dataset_meta.manual_seg == 1);
% dataset_meta = dataset_meta(segmented_id, :);
% dataset_bunny = dataset_meta(ismember(dataset_meta.stim_type, 'bunny'),:);
dataset_grat = dataset_meta(ismember(dataset_meta.stim_type, 'grating'),:);

%%

data = dataset_grat;
disp('analyzing grating datasets')

nset = size(data,1);
% for iset = 1:nset
% iset = 1 % session 002 220714 i1372 is all corrupted
iset = 2 % analyze session 003 220714 i1372

close all
clearvars -except dataset_grat iset nset data

iset, nset
dataset_now = data(iset,:);
arg_date = num2str(dataset_now.date)
arg_mouse = dataset_now.mouse
imouse = ['i', num2str(arg_mouse)];
area = 'V1'
sess_num = ['00', num2str(dataset_now.num)]

%%
% declare all global var for single dataset
global frame_rate range_base range_resp ...
    ncell ntrial trial_len_min ...
    nisi nori ori_list...
    id_ad id_noad id_isi2 id_isi3 id_ori ...

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

%% load

save_flag = 1; % toggle this to save/skip all .mat creation below

% get stim mat: input_behav for each session
xls_dir = fullfile(data_fn, imouse, arg_date)
cd(xls_dir)
xls_file = dir('*.xlsx');
data_now_meta = readtable(xls_file.name);

grat_id = find(contains(data_now_meta{:,1}, sess_num)); % use session number to find grating recording index
data_now_meta(grat_id,:)
frame_rate = data_now_meta.(5)(grat_id);

%%

for i = 1:length(grat_id)
    id = grat_id(i);
    time = data_now_meta.(8)(id);
    ImgFolder = data_now_meta.(1){id}(1:3); % '002_000_000' strip to '002'

    fName = fullfile(mworks_fn, ['data-' imouse '-' arg_date '-' num2str(time) '.mat']);
    temp = load(fName); % load behavior data "input", which clashes with built-in function
    input_behav_seq(i) = temp.input; clear temp
end

%%
areamousedate = [area '_' imouse '_' arg_date];
mapped_path = 'Z:\All_Staff\home\lan\Data\2P_images';

if dataset_now.cellpose_seg == 1
    result_folder = [mapped_path, '\mat_inter\', areamousedate, '_cellpose']; 
    disp('cellpose segm');
elseif dataset_now.manual_seg == 1
    result_folder = [mapped_path, '\mat_inter\', areamousedate]; 
    disp('manual segm');
end

if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

%% substitute npSub_tc w cellpose

cd(fullfile(tc_fn, [arg_date '_' imouse]))
cd([arg_date '_' imouse '_runs-', ImgFolder])

if dataset_now.cellpose_seg == 1
    tc = load([arg_date '_' imouse '_runs-', ImgFolder,'_TCs_cellpose.mat']);
elseif dataset_now.manual_seg == 1
    tc = load([arg_date '_' imouse '_runs-', ImgFolder,'_TCs_addfake.mat']);
end

disp('loaded cellpose or manual timecourse')
npSub_tc = tc.npSub_tc;
[nframe, ncell] = size(npSub_tc)
nframe_seq = [nframe];

%% params & indexing trials
% index by adapter contrast, target ori, isi

% % cut off corrupted frames: from frame 87139, count back 100 ms
% corrupt_buffer = frame_rate * 0.1 % 30 frame / s * 0.100 s = 3 frame
% frame_tg_onset = cell2mat(input_behav_seq.cTargetOn);
% ntrial_good = sum(frame_tg_onset < nframe - corrupt_buffer)
% 
% field_name = fieldnames(input_behav_seq);
% for i = 1 : length(field_name)
% %     field_content = getfield(input_behav_seq, field_name{i});
%     field_content = input_behav_seq.(field_name{i});
%     field_len = length(field_content);
% 
%     if field_len > ntrial_good
% %         input_behav_seq = setfield(input_behav_seq, field_name{i}, field_content{1:ntrial_good})
%         input_behav_seq.(field_name{i}) = field_content(1:ntrial_good);
%     end
% end

input_behav = input_behav_seq;
ntrial = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames
% ntrial = ntrial_good - 1;
% disp('chopped off corrupted data')

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

ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg); 
ori_seq(ori_seq == 180) = 0;
ori_seq(end) = [];
ori_list = unique(ori_seq); 
nori = length(ori_list); id_ori = cell(nori, 1);
for iori  = 1 : nori
    id_ori{iori} = find(ori_seq == ori_list(iori)); 
end

t = cellfun(@size,id_ori,'uni',false);
t = cell2mat(t(:,1));
nrep_stim = unique(t(:,2))

%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

cd(result_folder)

tc_align_ad = align_tc(frame_ad, npSub_tc);
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad);
dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad);

disp('extract only trace_by_trial for trials w isi=250, ad & tg ori=0')
trial_select = intersect(id_ori{find(ori_list==0)}, id_ad250);
trace_by_trial = dfof_align_ad(:, trial_select, :);
disp('stim_seq is just tg ori=0 for now')
stim_seq = ori_seq';
if save_flag; save trace_trial_stim.mat trace_by_trial stim_seq; end

%% set resp window
% find base window & resp window

t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); 
t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); 
t_tg = squeeze(nanmean(t(:,:), 1)); 

% figure
% range = 60;
% plot(t_ad(1:range), 'r'); hold on; 
% plot(t_tg(1:range), 'b'); 
% grid on; grid minor; 
% set(gcf, 'Position', get(0, 'Screensize')); 
% legend('ad align', 'targ align')
% if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end
% 
% % t_ad = [t_ad(100:end), t_ad(1:100)];
% % t_tg = [t_tg(100:end), t_tg(1:100)];
% % endpoint = length(t_ad);
% 
% figure
% range = trial_len_min; 
% plot(t_ad(1:range), 'r'); hold on; 
% plot(t_tg(1:range), 'b'); 
% % xline(endpoint - 100);
% grid on; grid minor; 
% set(gcf, 'Position', get(0, 'Screensize')); 
% legend('ad align', 'targ align')
% if save_flag; saveas(gcf, 'dfof align', 'jpg'); end

%% find calcium resp window

find_peak_bef = 15;
disp('assume: first peak comes before n frames. second peak comes after')
trace_start = t_ad(1:find_peak_bef);
[~, peak_id] = max(trace_start)
if peak_id < 6 % first peak should not be earlier than 6 frames
    disp('WARNING: strange trace or peak')
end
range_base = 1:3; range_resp = (peak_id-1):(peak_id+1);

figure;
% plot(t_ad(1:40));
plot(t_ad(1:end));
xline(range_base(end), 'k--')
xline(range_resp(1), 'k--')
xline(range_resp(end), 'k--')
grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['find_ca_latency_ca_window.jpg'])

%% trial-wise response and baseline

[dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align_ad, save_flag);
dfof_ad_trial = dfof_ad_trial(:,:,3); % only take isi=250
dfof_base_trial = dfof_base_trial(:,:,3);

[dfof_tg_trial, dfof_base2_trial] = dfof_resp_trialwise(dfof_align_tg, save_flag);
dfof_tg_trial = dfof_tg_trial(:,:,3);
dfof_base2_trial = dfof_base2_trial(:,:,3);

if save_flag; save resp_base_trialwise.mat dfof_ad_trial dfof_tg_trial...
        dfof_base_trial dfof_base2_trial; end

%% response to adapter & targets. get trace (bunny mode: isi=250 only)
% dfof_ad = ncell x nstim. dfof_tg = ncell x nstim

[dfof_ad, dfof_ad_sem, dfof_ad_std] = dfof_resp(dfof_align_ad, 'tg', 0); % tg mode aka separate diff stim images, but use adapter resp
[dfof_tg, dfof_tg_sem, dfof_tg_std] = dfof_resp(dfof_align_tg, 'tg', 0);
% dfof_ad = dfof_ad(:,:,3); dfof_ad_sem = dfof_ad_sem(:,:,3); dfof_ad_std = dfof_ad_std(:,:,3); % keep isi=250 only
% dfof_tg = dfof_tg(:,:,3); dfof_tg_sem = dfof_tg_sem(:,:,3); dfof_tg_std = dfof_tg_std(:,:,3);
if save_flag; save dfof.mat dfof_ad dfof_ad_sem dfof_ad_std dfof_tg dfof_tg_sem dfof_tg_std; end 

% trace = ncell x nori x nisi3 [noad 750 250]
[trace_avg, trace_sem] = trace_grand_avg(dfof_align_ad, 0);
size(trace_avg)
% trace_avg = squeeze(trace_avg(:,:,3,:)); trace_sem = squeeze(trace_sem(:,:,3,:));
if save_flag; save trace_aligned.mat trace_avg trace_sem; end

%% trial-wise response and baseline

[dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align_ad, save_flag);
% dfof_ad_trial = dfof_ad_trial(:,:,3);
% dfof_base_trial = dfof_base_trial(:,:,3);

[dfof_tg_trial, dfof_base2_trial] = dfof_resp_trialwise(dfof_align_tg, save_flag);
% dfof_tg_trial = dfof_tg_trial(:,:,3);
% dfof_base2_trial = dfof_base2_trial(:,:,3);

if save_flag; save resp_base_trialwise.mat dfof_ad_trial dfof_tg_trial...
        dfof_base_trial dfof_base2_trial; end

%
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