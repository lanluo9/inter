%% init
% this script takes in a MULTI session timecourse from one day, and aligns it to trial start and stimulus id
% output: trace_trial_stim.mat 
% code is inherited from adp_mag_pipeline_mix50.m and adp_mag_pipeline_grating.m

close all; clc; 
clear; clear global

root_path = 'C:\Users\GlickfeldLab\Documents\test\inter'; %'C:\Users\ll357\Documents\inter';
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dir_meta = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/adp_dataset_master.xlsx';
dataset_meta = readtable(dir_meta);

dataset_table = dataset_meta( ...
    logical(strcmp(dataset_meta.paradigm, 'grat_SF6') ...
    + strcmp(dataset_meta.paradigm, 'grat_SF_lindsey')), ...
    :);

dataset_table = dataset_table(~contains(dataset_table.note, 'bad'), :); % exclude bad

ndate = length(unique(dataset_table.date));
date_arr = unique(dataset_table.date);

%% find TC.mat for each date

for idate = 1:ndate

dataset_date = dataset_table(dataset_table.date == date_arr(idate), :);
clear global
close all

dataset_now = dataset_date(end, :); % get metadata from final sess of this date
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num{1}
area = dataset_now.area{1}
imouse = ['i', num2str(arg_mouse)];

dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
%     mkdir(dir_analysis)
    dir_analysis = strrep(dir_analysis, 'lan', 'lindsey')
    disp('this dir exist as a shortcut')
    cd(dir_analysis)
end

% % check if time course exists for final sess. only useful for multi day, multi sess data
if ~isempty(dir('*TC*.mat'))
    disp('timecourse exists, start generating trace_trial_stim.mat')
else
    disp('timecourse doesnt exist yet, skip to next set')
    continue
end

%%
% declare all global var for single dataset
global frame_rate range_base range_resp ...
    ncell ntrial trial_len_min ...
    nisi nori ori_list...
    id_ad id_noad id_isi2 id_isi3 id_ori

save_flag = 0; % toggle this to save/skip all .mat creation below

try % if data folder contains 2p imaging note.xls
    stim_protocol = 'grat_SF6_allen_nat8' % 'grat_SF6'
    xls_dir = fullfile(data_fn, imouse, arg_date)
    cd(xls_dir)
    xls_file = dir('*.xlsx');
    data_now_meta = readtable(xls_file.name);
    sess_id_arr = find(contains(data_now_meta{:,9}, stim_protocol));
    data_now_meta(sess_id_arr,:)
    frame_rate = data_now_meta.(5)(sess_id_arr(end));
catch % if data folder has been transferred to AWS, cannot read 2p imaging note.xls anymore
    sess_id_arr = str2num(dataset_now.num{1}(end));
    disp('only works for single session data. TODO: fix sess_id_arr and ImgFolder')
end

%% load input_behav & npSubTC for eash sess

df_flat = [];
nframe_seq = [];

for i = 1:length(sess_id_arr)
    ImgFolder = dataset_date.num{i};
    fName = fullfile(tc_fn, [arg_date, '_', imouse], [arg_date, '_', imouse, '_runs-', ImgFolder], ...
        [arg_date, '_', imouse, '_runs-', ImgFolder, '_input.mat']);
    temp = load(fName); % load behavior data "input", which clashes with built-in function
    try
        input_behav_seq(i) = temp.input; 
        frame_rate = temp.input.frameRateHz;
    catch
        input_behav_seq(i) = temp.behav_input; 
        frame_rate = temp.behav_input.frameRateHz;
    end
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
nframe, ncell, size(df_flat)

disp('loaded cellpose timecourse & visual input')
try
    assert(sum(df_flat(nframe,:) < 65535) == ncell)
catch
    plot(df_flat(nframe,:)) 
    disp('some cell at the final frame has fluo value 65535 -> corrupted')
end

areamousedate = [area '_' imouse '_' arg_date];
mapped_path = 'Z:\All_Staff\home\lan\Data\2P_images';
result_folder = [mapped_path, '\mat_inter\', areamousedate, segment_suffix];
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

if ~isempty(dir('*trace_trial_stim*.mat')) && ~isempty(dir('*resp_base_trialwise*.mat'))
    disp('TC align result exists, skip to next dataset:')
    continue
end

%% concat trial stim info for each session
% index by adapter contrast, target ori, isi

ntrial = 0;
adapter_id = [];
target_id = [];
frame_ad = [];
frame_ad_off = [];
frame_tg = [];

for isess = 1 : length(sess_id_arr)
    input_behav = input_behav_seq(isess);
    ntrial_sess = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames
    ntrial = ntrial + ntrial_sess;
    
    disp('only support same stim1-stim2 pairs')
    try % randStim1_doSameStim2 with bunnies6.mwel
        adapter_id_sess = cell2mat(input_behav.tstimOne);
%         target_id_sess = cell2mat(input_behav.tstimTwo);
    catch % grat_SF6 with twoStim.mwel
        SF_arr = sort(unique(cell2mat(input_behav.tStimOneGratingSpatialFreqCPD)));
        adapter_SF = cell2mat(input_behav.tStimOneGratingSpatialFreqCPD);
        adapter_id_sess = zeros(length(adapter_SF), 1);
        for iSF = 1 : length(adapter_SF)
            adapter_id_sess(iSF) = find(SF_arr==adapter_SF(iSF));
        end
        tmp = cell2mat(input_behav.tStimTwoGratingSpatialFreqCPD) == cell2mat(input_behav.tStimOneGratingSpatialFreqCPD);
        assert(sum(tmp) == length(tmp)) % stim 1 vs 2 have same SF
%         target_id_sess = adapter_id_sess;
    end

    adapter_id_sess = adapter_id_sess(1:ntrial_sess);
    try
        assert(size(adapter_id_sess, 1) > size(adapter_id_sess, 2)) % ensure its a column vector
    catch
        adapter_id_sess = adapter_id_sess';
    end
    adapter_id = [adapter_id; adapter_id_sess];

    target_id_sess = adapter_id_sess; % only support same stim1-stim2 pairs
    target_id = [target_id; target_id_sess];

    frame_ad_sess = double(cell2mat(input_behav.cStimOneOn)); 
    frame_ad_sess = frame_ad_sess(1:ntrial_sess);
%     frame_ad_off_sess = double(cell2mat(input_behav.cStimOneOn)); % TODO: fix bug!!!
%     frame_ad_off_sess = frame_ad_off_sess(1:ntrial_sess);
    frame_tg_sess = double(cell2mat(input_behav.cStimTwoOn)); 
    frame_tg_sess = frame_tg_sess(1:ntrial_sess);
    if isess > 1
        frame_ad_sess = frame_ad_sess + sum(nframe_seq(1:isess-1));
%         frame_ad_off_sess = frame_ad_off_sess + sum(nframe_seq(1:isess-1));
        frame_tg_sess = frame_tg_sess + sum(nframe_seq(1:isess-1));
    end
    frame_ad = [frame_ad, frame_ad_sess];
%     frame_ad_off = [frame_ad_off, frame_ad_off_sess];
    frame_tg = [frame_tg, frame_tg_sess];
end

adapter_list = unique(adapter_id); n_adapter = length(adapter_list);
target_list = unique(target_id); n_target = length(target_list);
id_noad = []; id_ad = 1:ntrial; 
% isi_seq = frame_tg - frame_ad_off; % isi_seq was never used, so bug in frame_ad_off_sess did not affect anything
trial_len_min = min(unique(diff(frame_ad)));

nisi = length(unique(frame_tg - frame_ad));
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

t = cellfun(@size,id_ori,'uni',false);
t = cell2mat(t(:,1));
nrep_stim = unique(t(:))
disp('ignore `1`')

%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

cd(result_folder)

npSub_tc = df_flat;
tc_align_ad = align_tc(frame_ad, npSub_tc);
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad);
dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad);

trace_by_trial = dfof_align_ad;
stim_seq = adapter_id';
if save_flag; save trace_trial_stim.mat trace_by_trial stim_seq; end

%% set resp window
% find base window & resp window

t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); 
t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); 
t_tg = squeeze(nanmean(t(:,:), 1)); 

figure
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
range = trial_len_min; 
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
trace_start = t_ad(1:find_peak_bef);
[~, peak_id] = max(trace_start)
if peak_id < 6 % first peak should not be earlier than 6 frames
    disp('WARNING: strange trace or peak')
end
range_base = 1:3
range_resp = (peak_id-1):(peak_id+2); % fast climb, slow fall
range_resp = range_resp - 1

figure;
plot(t_ad(1:40)); hold on;
plot(t_tg(1:40));
xline(range_base(end), 'k--')
xline(range_resp(1), 'k--')
xline(range_resp(end), 'k--')
grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'find_ca_latency_ca_window.jpg')

%% use resp window to cut trace_by_trial
% aka slicing dfof_align_ad & dfof_align_tg

close all

R1_cell_trial = mean(dfof_align_ad(:, :, range_resp), 3) ...
              - mean(dfof_align_ad(:, :, range_base), 3); % ncell x ntrial, avg over resp time win
R2_cell_trial = mean(dfof_align_tg(:, :, range_resp), 3) ...
              - mean(dfof_align_tg(:, :, range_base), 3);

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
dfof_ad = dfof_ad(:,:,3); dfof_ad_sem = dfof_ad_sem(:,:,3); dfof_ad_std = dfof_ad_std(:,:,3); % keep isi=250 only
dfof_tg = dfof_tg(:,:,3); dfof_tg_sem = dfof_tg_sem(:,:,3); dfof_tg_std = dfof_tg_std(:,:,3);
if save_flag; save dfof.mat dfof_ad dfof_ad_sem dfof_ad_std dfof_tg dfof_tg_sem dfof_tg_std; end 

% trace = ncell x nori x nisi3 [noad 750 250]
[trace_avg, trace_sem] = trace_grand_avg(dfof_align_ad, 0);
trace_avg = squeeze(trace_avg(:,:,3,:)); % trace_sem = squeeze(trace_sem(:,:,3,:));
if save_flag; save trace_aligned.mat trace_avg trace_sem; end

%% trial-wise response and baseline

[dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align_ad, save_flag);
dfof_ad_trial = dfof_ad_trial(:,:,3);
dfof_base_trial = dfof_base_trial(:,:,3);

[dfof_tg_trial, dfof_base2_trial] = dfof_resp_trialwise(dfof_align_tg, save_flag);
dfof_tg_trial = dfof_tg_trial(:,:,3);
dfof_base2_trial = dfof_base2_trial(:,:,3);

if save_flag; save resp_base_trialwise.mat dfof_ad_trial dfof_tg_trial...
        dfof_base_trial dfof_base2_trial; end

% %% discard trials by pupil or run speed -> trial_filter_by_pupil_or_speed.m
% open trial_filter_by_pupil_or_speed.m
% % TODO: when we integrate trial filter (bool) into dataframe, need to cut
% % off final trial of each sess just like above:
% % ntrial_sess = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames

%% find visually driven cells -> vis_driven.ipynb

end