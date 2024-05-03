%% init
% this script takes in a single session timecourse from one day, 
% and aligns it to trial start and stimulus id

% NOTE: there can be multiple single sessions (typically diff depth of 
% same mouse & area) from one single day, so outer loop is over date_sess, not date alone

% output: trace_trial_stim.mat 
% code is inherited from adp_mag_pipeline_mix50.m and adp_mag_pipeline_grating.m

close all; clc; 
clear; clear global

root_path = 'C:\Users\ll357\Documents\inter';
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dir_meta = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/adp_dataset_master.xlsx';
dataset_meta = readtable(dir_meta);

stim_type = 'grating' % grat_8ori_3isi
dataset_table = dataset_meta(strcmp(dataset_meta.paradigm, stim_type), :);

dataset_table = dataset_table(strcmp(dataset_table.gcamp, '6s'), :);
seg_bool = dataset_table.manual_seg | dataset_table.cellpose_seg; % exclude not-segmented data
dataset_table = dataset_table(seg_bool, :);

% area_bool = logical(strcmp(dataset_table.area, 'LM') + strcmp(dataset_table.area, 'V1'));
% area_bool = logical(strcmp(dataset_table.area, 'LI'));
% dataset_table = dataset_table(area_bool, :);
sum(strcmp(dataset_table.area, 'V1'))
sum(strcmp(dataset_table.area, 'LM'))
sum(strcmp(dataset_table.area, 'LI'))

% dataset_table = dataset_table(dataset_table.date == 210120, :)

nset = size(dataset_table, 1);

%% find TC.mat

for iset = 1:nset

% iset = 1
clear global; close all
iset, nset

dataset_now = dataset_table(iset,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num{1}
area = dataset_now.area{1}

imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end

% check if time course exists
disp('wait for timecourse')
while isempty(dir('*TC*.mat'))
    pause(60) % wait for TC extract to finish
end
disp('timecourse exists, start generating trace_trial_stim.mat')

% if ~isempty(dir('*TC*.mat'))
%     disp('timecourse exists, start generating trace_trial_stim.mat')
% else
%     disp('timecourse doesnt exist yet, waiting for it')
%     pause(60) % wait for TC extract to finish
% end

%%
% declare all global var for single dataset
global frame_rate range_base range_resp ...
    ncell ntrial trial_len_min ...
    nisi nori ori_list...
    id_ad id_noad id_isi2 id_isi3 id_ori

save_flag = 0; % toggle this to save/skip all .mat creation below

% try % if data folder contains 2p imaging note.xls
%     stim_protocol = 'grat_SF6'
%     xls_dir = fullfile(data_fn, imouse, arg_date)
%     cd(xls_dir)
%     xls_file = dir('*.xlsx');
%     data_now_meta = readtable(xls_file.name);
%     sess_id_arr = find(contains(data_now_meta{:,9}, stim_protocol));
%     data_now_meta(sess_id_arr,:)
%     frame_rate = data_now_meta.(5)(sess_id_arr(end));
% catch % if data folder has been transferred to AWS, cannot read 2p imaging note.xls anymore
    sess_id_arr = str2num(dataset_now.num{1}(end));
    disp('only works for single session data. TODO: fix sess_id_arr and ImgFolder')
% end

%% load input_behav & npSubTC for eash sess

df_flat = [];
nframe_seq = [];
clear input_behav_seq

for i = 1:length(sess_id_arr)
    ImgFolder = dataset_now.num{i};
    fName = fullfile(tc_fn, [arg_date, '_', imouse], [arg_date, '_', imouse, '_runs-', ImgFolder], ...
        [arg_date, '_', imouse, '_runs-', ImgFolder, '_input.mat']);
    temp = load(fName); % load behavior data "input", which clashes with built-in function
    try
        input_behav_seq(i) = temp.input; 
    catch
        input_behav_seq(i) = temp.behav_input;
    end

    frame_rate = input_behav_seq(i).frameRateHz;
    clear temp

    cd(fullfile(tc_fn, [arg_date '_' imouse]))
    cd([arg_date '_' imouse '_runs-', ImgFolder])
    segment_suffix = ''; % default suffix is empty -> using manual segment
    try
        tc = load([arg_date '_' imouse '_runs-', ImgFolder,'_TCs_addfake.mat']); % try manual TC first
    catch
        tc = load([arg_date '_' imouse '_runs-', ImgFolder,'_TCs_cellpose.mat']); % if not, use cellpose TC
        segment_suffix = '_cellpose'; % add cellpose suffix for mat_inter subdir, if segmented by cellpose
    end
    df_flat = [df_flat; tc.npSub_tc];

    [nframe, ncell] = size(tc.npSub_tc);
    nframe_seq = [nframe_seq, nframe];
end
nframe, ncell

disp('loaded timecourse & visual input')
try
    assert(sum(df_flat(nframe,:) < 65535) == ncell)
catch
    plot(df_flat(nframe,:)) 
    disp('some cell at the final frame has fluo value 65535 -> corrupted')
end

area_mouse_date_sess = [area '_' imouse '_' arg_date '_' arg_ImgFolder];
mapped_path = 'Z:\All_Staff\home\lan\Data\2P_images';
result_folder = [mapped_path, '\mat_inter\', area_mouse_date_sess, segment_suffix];
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

%% find trial stim info for each session
% index by adapter contrast, target ori, isi

input_behav = input_behav_seq;
ntrial = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames

contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast); 
id_noad = find(contrast_ad == 0); id_ad = find(contrast_ad == 1); 
id_noad(id_noad > ntrial) = []; id_ad(id_ad > ntrial) = []; 

frame_ad = double(cell2mat(input_behav.cStimOn)); frame_ad_off = double(cell2mat(input_behav.cStimOff));
frame_tg = celleqel2mat_padded(input_behav.cTargetOn); frame_tg = double(frame_tg);
isi_seq = frame_tg - frame_ad_off; 

% nisi = length(unique(frame_tg - frame_ad)); % NOTE: fixed bug 230704, did
% not rerun pipeline bc it doens't seem to affect anything

id_750 = find(isi_seq > mean(isi_seq)); id_250 = find(isi_seq < mean(isi_seq)); 
id_750(id_750 > ntrial) = []; id_250(id_250 > ntrial) = []; 
id_ad750 = intersect(id_ad, id_750); id_ad250 = intersect(id_ad, id_250);
id_isi2 = {id_ad750, id_ad250}; 
id_isi3 = {id_noad, id_ad750, id_ad250};
nisi = length(id_isi3); % redefine n_isi
trial_len_min = min(unique(diff(frame_ad)));

paradigm_ms.stim1_ms = input_behav.stimOnTimeMs;
paradigm_ms.stim2_ms = input_behav.targetOnTimeMs;
paradigm_ms.max_isi_ms = max(isi_seq) / frame_rate * 1000; % fixed frame-to-ms bug on 20230703
paradigm_ms.iti_ms = input_behav.itiTimeMs;

ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg); 
ori_seq(ori_seq == 180) = 0;
ori_seq(end) = [];
ori_list = unique(ori_seq); 
nori = length(ori_list); 
id_ori = cell(nori, 1);
for iori  = 1 : nori
    id_ori{iori} = find(ori_seq == ori_list(iori)); 
end

t = cellfun(@size,id_ori,'uni',false);
t = cell2mat(t(:,1));
nrep_stim = unique(t(:,2))

%% dfof aligned
% align tc by adapter or targ onset. normalize by "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

cd(result_folder)

% if ~isempty(dir('*trace_trial_stim*.mat')) && ~isempty(dir('*resp_base_trialwise*.mat'))
%     disp('TC align result exists, skip to next dataset:')
%     continue
% end

npSub_tc = df_flat; % nframe x ncell
tc_align_ad = align_tc(frame_ad, npSub_tc); % ncell x ntrial x nframe_trial
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad, paradigm_ms); % same as above but df/f
dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad, paradigm_ms);

trace_by_trial = dfof_align_ad;
stim_ori = ori_seq'; % stim as col
isi_nframe = isi_seq'; % ISI as number of frames in each trial
adapter_contrast = contrast_ad'; % contrast of adapter (R1)
if save_flag; save trace_trial_stim.mat trace_by_trial ...
        stim_ori isi_nframe adapter_contrast; end

% %% san check
% % what can generate possibly fake adp when adapter vs target are orthogonal?  
% % check if bin=90 tuning curve is real: plot timecourse for stim2=90, noad vs 250
% 
% size(dfof_align_ad) % ncell x ntrial x nframe
% 
% isi_nframe = isi_nframe(1:length(stim_ori));
% adapter_contrast = adapter_contrast(1:length(stim_ori));
% stim_id_now = 90-22.5;
% trial_id_noad_1 = (stim_ori==stim_id_now) & (isi_nframe<10) & (adapter_contrast==0); % trials without adapter but fake isi=250, stim2 ori=90
% trial_id_noad_2 = (stim_ori==stim_id_now) & (isi_nframe>10) & (adapter_contrast==0); % trials without adapter but fake isi=750, stim2 ori=90
% trial_id_ad_250 = (stim_ori==stim_id_now) & (isi_nframe<10) & (adapter_contrast==1); % trials with isi=250 adapter, stim2 ori=90
% trial_id_ad_750 = (stim_ori==stim_id_now) & (isi_nframe>10) & (adapter_contrast==1); % trials with isi=250 adapter, stim2 ori=90
% sum(trial_id_noad_1)
% sum(trial_id_noad_2)
% sum(trial_id_ad_250)
% sum(trial_id_ad_750)
% 
% file = 'C:\Users\ll357\Documents\inter\results\tuning curve bias san check\vis_driven_ttest_bonferroni.mat';
% tmp = load(file);
% cell_id_filter = logical(tmp.vis_driven');
% cell_id_filter_new = logical(tmp.vis_driven');
% clear tmp
% 
% % file = 'C:\Users\ll357\Documents\inter\results\tuning curve bias san check\vis_cell_bool.mat'; % old vis cells, from anova/kruskal AND evoked thresh
% % tmp = load(file);
% % cell_id_filter = logical(tmp.vis_cell_bool');
% % cell_id_filter_old = logical(tmp.vis_cell_bool');
% % clear tmp
% 
% % cell_id_filter = ~cell_id_filter; % only plot not-vis-driven cells
% cell_id_filter = logical(ones(size(cell_id_filter))); % turn off cell filter
% % cell_id_filter = ones(size(cell_id_filter), 'like', cell_id_filter); % turn off cell filter
% 
% 
% trace_noad_1 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_noad_1, :), 1), 2));
% trace_noad_2 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_noad_2, :), 1), 2));
% trace_ad_250 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_ad_250, :), 1), 2));
% 
% % trace_noad_1 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_1, :), 1), 2));
% % trace_noad_2 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_2, :), 1), 2));
% % trace_ad_250 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_ad_250, :), 1), 2));
% 
% figure;
% plot(trace_noad_1, 'b')
% hold on
% plot(trace_noad_2, 'g')
% plot(trace_ad_250, 'r')
% legend('noad early', 'noad late', 'ad250');
% xlim([0, 60]);
% 
% 
% % vis_driven = cell_id_filter;
% % for icell = 1 : 10 %size(cell_id_filter, 1)
% % %     if vis_driven(icell) == 1
% % %         continue
% % %     end
% %     cell_id_filter = zeros(size(cell_id_filter));
% %     cell_id_filter(icell) = 1;
% %     cell_id_filter = logical(cell_id_filter); % only use a single cell
% %     
% %     trace_noad_1 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_1, :), 1), 2));
% %     trace_noad_2 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_2, :), 1), 2));
% %     trace_ad_250 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_ad_250, :), 1), 2));
% %     
% %     figure;
% %     plot(trace_noad_1, 'b')
% %     hold on
% %     plot(trace_noad_2, 'g')
% %     plot(trace_ad_250, 'r')
% %     legend('noad early', 'noad late', 'ad250');
% %     xlim([0, 60]);
% %     title(num2str(vis_driven(icell)))
% %     pause(2)
% % end

%% set resp window
% find base window & resp window

t = squeeze(nanmean(dfof_align_ad, 1));
t_ad = squeeze(nanmean(t, 1)); 
t = squeeze(nanmean(dfof_align_tg, 1));
t_tg = squeeze(nanmean(t, 1)); 

figure;
range = 60;
plot(t_ad(1:range), 'r'); hold on; 
plot(t_tg(1:range), 'b'); 
grid on; grid minor; 
set(gcf, 'Position', get(0, 'Screensize')); 
legend('ad align', 'targ align')
if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end

% t_ad = [t_ad(100:end), t_ad(1:100)];
% t_tg = [t_tg(100:end), t_tg(1:100)];
% endpoint = length(t_ad);

figure
range = trial_len_min-1; 
plot(t_ad(1:range), 'r'); hold on; 
plot(t_tg(1:range), 'b'); 
% xline(endpoint - 100);
grid on; grid minor; 
set(gcf, 'Position', get(0, 'Screensize')); 
legend('ad align', 'targ align')
if save_flag; saveas(gcf, 'dfof align', 'jpg'); end

%% find calcium resp window

find_peak_bef = 15;
disp('assume: first peak comes before n frames. second peak comes after')
trace_start = t_ad(1:find_peak_bef) + t_tg(1:find_peak_bef);
[~, peak_id] = max(trace_start)
if peak_id < 6 % first peak should not be earlier than 6 frames
    disp('WARNING: strange trace or peak')
end
range_base = 1:5; range_resp = (peak_id-2):(peak_id+2);

figure;
plot(t_ad(1:40)); hold on;
plot(t_tg(1:40));
xline(range_base(end), 'k--')
xline(range_resp(1), 'k--')
xline(range_resp(end), 'k--')
grid on; grid minor
% set(gcf, 'Position', get(0, 'Screensize'));
if save_flag; saveas(gcf, 'find_ca_latency_ca_window.jpg'); end

%% use resp window to cut trace_by_trial
% aka slicing dfof_align_ad & dfof_align_tg

close all

R1_cell_trial = mean(dfof_align_ad(:, :, range_resp), 3) ...
              - mean(dfof_align_ad(:, :, range_base), 3); % ncell x ntrial, avg over resp time win
R2_cell_trial = mean(dfof_align_tg(:, :, range_resp), 3) ...
              - mean(dfof_align_tg(:, :, range_base), 3); % separate 250 vs 750 isi by saved var isi_nframe

figure
subplot(211)
imshow(R1_cell_trial, 'InitialMagnification', 800);
colorbar;
subplot(212)
imshow(R2_cell_trial, 'InitialMagnification', 800);
colorbar;
set(gcf, 'Position', get(0, 'Screensize'));

if save_flag; save('trace_trial_stim.mat', 'R1_cell_trial', 'R2_cell_trial', '-append'); end

%% response to adapter & targets. get trace (bunny mode: isi=250 only)
% dfof_ad = ncell x nstim. dfof_tg = ncell x nstim

close all
[dfof_ad, dfof_ad_sem, dfof_ad_std] = dfof_resp(dfof_align_ad, 'tg', 0); % tg mode aka separate diff stim images, but use adapter resp
[dfof_tg, dfof_tg_sem, dfof_tg_std] = dfof_resp(dfof_align_tg, 'tg', 0);
% dfof_ad = dfof_ad(:,:,3); dfof_ad_sem = dfof_ad_sem(:,:,3); dfof_ad_std = dfof_ad_std(:,:,3); % keep isi=250 only
% dfof_tg = dfof_tg(:,:,3); dfof_tg_sem = dfof_tg_sem(:,:,3); dfof_tg_std = dfof_tg_std(:,:,3);
if save_flag; save dfof.mat dfof_ad dfof_ad_sem dfof_ad_std dfof_tg dfof_tg_sem dfof_tg_std; end 

% trace = ncell x nori x nisi3 [noad 750 250]
[trace_avg, trace_sem] = trace_grand_avg(dfof_align_ad, 0);
% trace_avg = squeeze(trace_avg(:,:,3,:)); % trace_sem = squeeze(trace_sem(:,:,3,:));
if save_flag; save trace_aligned.mat trace_avg trace_sem; end

% %% san check w only cells preferring stim2=90
% 
% dfof_tg_noad = dfof_tg(:, :, 1); % no adapter, ncell x nstim
% [argvalue, argmax] = max(dfof_tg_noad, [], 2);
% cell_pref_90 = logical(argmax == 5); % 0, 22.5 .. 90 -> stim id = 5 is 90. NOTE: 1-based indexing!!!
% % cell_pref_almost_0 = logical((argmax == 8) + (argmax == 1) + (argmax == 0)); % 0, 22.5 .. 90 -> stim id = 0 is 0 deg
% 
% cell_id_filter = cell_pref_90;
% sum(cell_id_filter)
% dfof_tg_noad_cell_pref = nanmean(dfof_tg_noad(cell_id_filter, :), 1);
% figure; plot(dfof_tg_noad_cell_pref)
% 
% trace_noad_1 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_noad_1, :), 1), 2));
% trace_noad_2 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_noad_2, :), 1), 2));
% trace_ad_250 = squeeze(nanmean(nanmean(dfof_align_ad(cell_id_filter, trial_id_ad_250, :), 1), 2));
% 
% % trace_noad_1 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_1, :), 1), 2));
% % trace_noad_2 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_noad_2, :), 1), 2));
% % trace_ad_250 = squeeze(nanmedian(nanmedian(dfof_align_ad(cell_id_filter, trial_id_ad_250, :), 1), 2));
% 
% figure;
% plot(trace_noad_1, 'b')
% hold on
% plot(trace_noad_2, 'g')
% plot(trace_ad_250, 'r')
% legend('noad early', 'noad late', 'ad250');
% xlim([0, 60]);


%% trial-wise response and baseline

[dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align_ad, save_flag);
% dfof_ad_trial = dfof_ad_trial(:,:,3);
% dfof_base_trial = dfof_base_trial(:,:,3);

[dfof_tg_trial, dfof_base2_trial] = dfof_resp_trialwise(dfof_align_tg, save_flag);
% dfof_tg_trial = dfof_tg_trial(:,:,3);
% dfof_base2_trial = dfof_base2_trial(:,:,3);

if save_flag; save resp_base_trialwise.mat dfof_ad_trial dfof_tg_trial...
        dfof_base_trial dfof_base2_trial; end

%% trial-wise response and baseline for Jeff population vector decoder

[ppResp] = dfof_resp_trialwise_jeff(dfof_align_tg, save_flag);

% %% discard trials by pupil or run speed -> trial_filter_by_pupil_or_speed.m
% open trial_filter_by_pupil_or_speed.m
% TODO: when we integrate trial filter (bool) into dataframe, need to cut
% off final trial of each sess just like above:
% ntrial_sess = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames

%% fit von mises tuning curve
% fit_param under conditions = ncell x nparam x nisi [noad vs ad750 vs ad250]
% ori_pref under conditions = ncell x nisi [noad vs ad750 vs ad250]
% save tuning curve for each cell & isi

[fit_param, ori_pref, tuning_curve_cell_isi] = fit_tuning(dfof_tg, save_flag);
[R_sq, ori_fit] = fit_tuning_jeff(dfof_tg, save_flag);

% %% fit tuning curve with half trials
% 
% % save_flag = 1
% fit_tuning_half_trials(dfof_align_tg, save_flag);
% save_flag = 0

%% well-fit cells
% cells whose noad-tg 90% bootstraps are within 22.5 deg of all-trials-included fit

bootstrap_file = fullfile(result_folder, 'fit_bootstrap_90perc.mat');
if exist(bootstrap_file, 'file')
    tmp = load(bootstrap_file, 'ori_perc');
    theta_90 = tmp.ori_perc'; % shape = [1 x ncell]. 90 percentile distance from avg fitted ori_pref
    
    % tmp = load(bootstrap_file, 'ori_pref_runs');
    % if size(tmp.ori_pref_runs, 2) == 1000
    %     disp('already done 1k bootstraps for well_fit, skip')
    % end
else
    save_flag = 1;
    nrun = 1000; 
    well_fit_cell = well_fit_cell_criteria(dfof_align_tg, nrun, save_flag); 

    tmp = load(bootstrap_file, 'ori_perc');
    theta_90 = tmp.ori_perc'; % shape = [1 x ncell]. 90 percentile distance from avg fitted ori_pref
end

% % relax well_fit_cell_criteria to a lower percentile for LI
% bootstrap_file = fullfile(result_folder, 'fit_bootstrap_relax.mat');
% if exist(bootstrap_file, 'file')
%     continue
% else
%     save_flag = 1;
%     nrun = 1000;
%     percentile_threshold = 0.7;
%     well_fit_cell = well_fit_cell_criteria_relax(percentile_threshold, nrun, save_flag);
%     disp('well_fit_cell percent:')
%     sum(well_fit_cell) / length(well_fit_cell) * 100
% end

% % % validation that well_fit cells' tuning is comparable no matter how you
% % % calculate their tuning: by avg or med of ori_pref_runs, or by ori_pref
% % % generated from fit_tuning func
% % tmp = load(bootstrap_file);
% % well_fit_bool = tmp.well_fit_cell;
% % 
% % ori_pref_noad = ori_pref(:, 1); % first col is no-adapter ori pref
% % ori_pref_noad_boot_avg = mean(tmp.ori_pref_runs, 2); % avg across boots
% % ori_pref_noad_boot_med = median(tmp.ori_pref_runs, 2);
% % 
% % ori_pref_noad = ori_pref_noad(well_fit_bool);
% % ori_pref_noad_boot_avg = ori_pref_noad_boot_avg(well_fit_bool);
% % ori_pref_noad_boot_med = ori_pref_noad_boot_med(well_fit_bool);
% % 
% % plot(ori_pref_noad_boot_avg)
% % hold on
% % plot(ori_pref_noad_boot_med)
% % plot(ori_pref_noad)
% % legend('boot avg', 'boot med', 'fit')

%% well max cells (limit 10)

well_max_file = fullfile(result_folder, 'well_max_10cell.mat');
try
    tmp = load(well_max_file, 'well_max');
    well_max = tmp.well_max; % shape = [1 x ncell]
catch
    disp('not enough vis cells, set well max to all false')
    well_max = zeros([1, ncell]);
end

%% find visually driven cells -> vis_driven.ipynb
% read pickle data

vis_file = fullfile(result_folder, 'vis_driven_ttest_bonferroni_jeff.mat');
tmp = load(vis_file);
vis_bool = tmp.vis_driven'; % column vector

%% save data for jeff population vector decoder, (un)masked with NaN

% size(ori_fit) % 181 x ncell
% size(R_sq) % ncell x 1
% size(theta_90) % 1 x ncell

% % ori_fit(:, ~vis_bool) = NaN;
% % R_sq(~vis_bool) = NaN;
% % theta_90(~vis_bool') = NaN;
% % well_max(~vis_bool') = NaN;

save_flag = 1;

if save_flag; save pop_vec_decoder_jeff_10wellmax_notnorm.mat ...
    ppResp ori_fit R_sq theta_90 well_max; end

% %% normalize pop vector
% 
% % % slice everything with well_max_10
% well_max_bool = logical(well_max);
% ppResp = cellfun(@(x) x(well_max_bool, :), ...
%                 ppResp, ...
%                 'UniformOutput',false);
% % ori_fit = ori_fit(:, well_max_bool);
% % R_sq = R_sq(well_max_bool);
% theta_90 = theta_90(well_max_bool);
% well_max = well_max(well_max_bool);
% 
% % % % for cell, min-max ori_fit to 0-1 range -> turned off, this seem to cause pv decoder perf = 0
% % ncell_sess = size(ori_fit, 2);
% % for icell = 1 : ncell_sess
% %     tuning_icell = ori_fit(:, icell);
% %     tuning_icell = (tuning_icell - min(tuning_icell)) / (max(tuning_icell) - min(tuning_icell));
% %     ori_fit(:, icell) = tuning_icell;
% % end
% 
% % % min-max doenst work & only ends up with AUROC full of zeros
% % % switch to re-fit ori_fit using normed dfof_tg
% dfof_tg = dfof_tg(well_max_bool, :, :);
% nori_sess = size(dfof_tg, 2);
% nisi_sess = size(dfof_tg, 3);
% for iisi = 1 : nisi_sess 
%     for iori = 1 : nori_sess
%         pop_vec_itrial = dfof_tg(:, iori, iisi); % pop vector in dfof_tg of isi and ori
%         pop_vec_normed = pop_vec_itrial / norm(pop_vec_itrial);
%         dfof_tg(:, iori, iisi) = pop_vec_normed;
%     end
% end
% [R_sq, ori_fit] = fit_tuning_jeff(dfof_tg, save_flag);
% 
% 
% % % for iisi x iori of ppResp, for itrial, get [ncell x 1] vector. normalize vector to len=1
% nisi_sess = size(ppResp, 1);
% nori_sess = size(ppResp, 2);
% for iisi = 1 : nisi_sess 
%     for iori = 1 : nori_sess
%         ntrial_sess = size(ppResp{iisi, iori}, 2);
%         for itrial = 1 : ntrial_sess
%             pop_vec_itrial = ppResp{iisi, iori}(:, itrial); % pop vector of single trial
%             pop_vec_normed = pop_vec_itrial / norm(pop_vec_itrial);
%             ppResp{iisi, iori}(:, itrial) = pop_vec_normed;
%         end
%     end
% end
% 
% if save_flag; save pop_vec_decoder_jeff_10wellmax_vecnorm.mat ...
%     ppResp ori_fit R_sq theta_90 well_max; end


%% (un)comment loop over sess

end
