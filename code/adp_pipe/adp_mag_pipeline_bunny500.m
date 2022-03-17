%% init

close all; clc; clear
clear id_ad id_noad id_isi2 id_isi3 id_ori root_path
clear frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list
global id_ad id_noad id_isi2 id_isi3 id_ori % declare all global var for single dataset
global frame_rate range_base range_resp ncell ntrial trial_len_min nisi nori ori_list
global root_path

root_path = 'C:\Users\ll357\Documents\inter';
cd([root_path, '\mat'])

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

% caiman bunny500/top gcamp6s
dataset_list = struct;
dataset_list.mouse = [1369];
dataset_list.date = [220310];
dataset_list.area = {'V1'};
stim_protocol = 'bunny'

%% load [xls, timecourse, stim]

save_flag = 1; % toggle this to save/skip all .mat creation below
iset = 1 % this assumes registered multisession is recorded on same day & mouse
date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse]
area = dataset_list.area{1,iset}

% get stim mat: input_behav for each session
xls_dir = fullfile(data_fn, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx');
data_now_meta = readtable(xls_file.name);
frame_rate = data_now_meta.(5)(end);
% bunny500_id = find(contains(data_now_meta.StimulusSet,stim_protocol));
% bunny500_id = find(contains(data_now_meta.Var9,stim_protocol));
bunny500_id = find(contains(data_now_meta.adjustIsi_500_StimDuration_200,stim_protocol));
data_now_meta(bunny500_id,:)

%%

clear input_behav_seq
for i = 1:length(bunny500_id)
    id = bunny500_id(i);
    time = data_now_meta.(8)(id);
    ImgFolder = data_now_meta.(1){id}(1:3);

    fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' num2str(time) '.mat']);
    temp = load(fName); % load behavior data "input"
    input_behav_seq(i) = temp.input; clear temp
end

sess_flag = ''
if sess_flag == 1
    input_behav_seq = input_behav_seq(1); sess = '_002';
elseif sess_flag == 3
    input_behav_seq = input_behav_seq(3); sess = '_004';
elseif contains(sess_flag, 'A')
    sess = '_sideA';
elseif contains(sess_flag, 'B')
    sess = '_sideB';
else
    sess = '';
end

areamousedate = [area '_' imouse '_' date sess];
mapped_path = 'Z:\All_Staff\home\lan\Data\2P_images';
% result_folder = [mapped_path, '\mat_inter\', areamousedate]; disp('manual segm');
result_folder = [mapped_path, '\mat_inter\', areamousedate, '_caiman']; disp('caiman segm');
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

%% substitute npSub_tc w caiman

% cd(fullfile(tc_fn, [date '_' imouse]))
% cd('220225_i1350_runs-002')
% tc = load('220225_i1350_runs-002_TCs_addfake.mat');
% df_flat = tc.npSub_tc;
% [nframe, ncell] = size(df_flat)

tc_file = fullfile(tc_fn, [date '_' imouse], ['caiman_activity_' imouse '_' date, '_multisess.mat']);
df = load(tc_file);
df_pile = (df.df)';

t = cellfun(@size,df_pile,'uni',false); 
ncell = size(t,2)
t = cell2mat(t(:,1));
nframe_seq = t(:,2);

df_flat = zeros(sum(nframe_seq), ncell); % [nframe_sum, ncell]
for icell = 1:ncell
    df_flat(:,icell) = horzcat(df_pile{:, icell})';
end

% if sess_flag == 1
%     frame_range = 1:70000;
% elseif sess_flag == 3
%     frame_range = 140000:210000;
% else
%     frame_range = 1:210000;
% end
% df_flat = df_flat(frame_range, :);
size(df_flat)

%% concat trial stim info for each session
% index by adapter contrast, target ori, isi

ntrial = 0;
adapter_id = [];
target_id = [];
frame_ad = [];
frame_ad_off = [];
frame_tg = [];

for i = 1 : length(bunny500_id) % comment out i>1 for single sess
    input_behav = input_behav_seq(i);
    ntrial_sess = input_behav.trialSinceReset - 1; % final trial discarded bc too few frames
    ntrial = ntrial + ntrial_sess;
    
    adapter_id_sess = cell2mat(input_behav.tstimOne); 
    adapter_id_sess = adapter_id_sess(1:ntrial_sess);
    adapter_id = [adapter_id, adapter_id_sess];
    target_id_sess = cell2mat(input_behav.tstimTwo); 
    target_id_sess = target_id_sess(1:ntrial_sess);
    target_id = [target_id, target_id_sess];
    
    assert(input_behav.doRandStimOne == 1 & input_behav.doSameStims == 1) % verify randStim1_doSameStim2 protocol
    assert(sum(adapter_id == target_id) == length(target_id))

    frame_ad_sess = double(cell2mat(input_behav.cStimOneOn)); 
    frame_ad_sess = frame_ad_sess(1:ntrial_sess);
    frame_ad_off_sess = double(cell2mat(input_behav.cStimOneOn)); 
    frame_ad_off_sess = frame_ad_off_sess(1:ntrial_sess);
    frame_tg_sess = double(celleqel2mat_padded(input_behav.cStimTwoOn)); 
    frame_tg_sess = frame_tg_sess(1:ntrial_sess);
    if i>1
        frame_ad_sess = frame_ad_sess + sum(nframe_seq(1:i-1));
        frame_ad_off_sess = frame_ad_off_sess + sum(nframe_seq(1:i-1));
        frame_tg_sess = frame_tg_sess + sum(nframe_seq(1:i-1));
    end
    frame_ad = [frame_ad, frame_ad_sess];
    frame_ad_off = [frame_ad_off, frame_ad_off_sess];
    frame_tg = [frame_tg, frame_tg_sess];
end

adapter_list = unique(adapter_id); n_adapter = length(adapter_list);
target_list = unique(target_id); n_target = length(target_list);
id_noad = []; id_ad = 1:ntrial; 
isi_seq = frame_tg - frame_ad_off; 
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
nrep_stim = unique(t(:,2)) 
% bunny500 3sess = 3/4/5 rep of each img
% bunnytop 3sess = 48-50 rep of each img, each sess = 16-17 rep

% %% visliz F over time: photobleaching / dying cell
% % clear all global % sbxread might need clear global var
% 
% cd('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002')
% names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002\*_000_000.mat');
% imgMatFile = {names.name}
% load(imgMatFile{1});
% nframes = info.config.frames;
% names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002\*_000_000.sbx');
% sbxFile = {names.name}
% data_temp = sbxread(sbxFile{1}(1:end-4), 0, nframes);
% 
% cd('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\003')
% names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\003\*_000_000.sbx');
% sbxFile = {names.name}
% data_temp2 = sbxread(sbxFile{1}(1:end-4), 0, nframes);
% 
% data = cat(4, data_temp, data_temp2);
% data = squeeze(data);
% data_avg = squeeze(mean(mean(data,2),1));
% 
% data_smooth = movmean(data_avg, 500);
% plot(data_smooth)
% 
% cd('Z:\All_Staff\home\lan\Analysis\2P\220225_i1350')
% save overall_fluorescence.mat data_avg
% saveas(gcf, 'overall_fluorescence', 'jpg')
% close

%% dfof aligned
% align tc by adapter or targ onset. normalize by 1-sec "trial baseline" to get dfof
% always use frame_ad as the end point of trial-specific baseline

npSub_tc = df_flat;
tc_align_ad = align_tc(frame_ad, npSub_tc);
tc_align_tg = align_tc(frame_tg, npSub_tc);
dfof_align_ad = tc_align_ad; % did not do trial-specific baselining
dfof_align_tg = tc_align_tg; % which should have been dfof_align = (tc - base) / base
% dfof_align_ad = dfof_by_trial_base(tc_align_ad, npSub_tc, frame_ad);
% dfof_align_tg = dfof_by_trial_base(tc_align_tg, npSub_tc, frame_ad);

trace_by_trial = dfof_align_ad;
stim_seq = adapter_id';
if save_flag; save trace_trial_stim.mat trace_by_trial stim_seq; end

%% set resp window
% find base window & resp window

t = squeeze(nanmean(squeeze(dfof_align_ad(:,:,:)), 1)); t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(dfof_align_tg(:,:,:)), 1)); t_tg = squeeze(nanmean(t(:,:), 1)); 

figure
range = 50;
plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize')); legend('ad align', 'targ align')
if save_flag; saveas(gcf, 'dfof align zoomin', 'jpg'); end

t_ad = [t_ad(100:end), t_ad(1:100)];
t_tg = [t_tg(100:end), t_tg(1:100)];
endpoint = length(t_ad)

figure
range = trial_len_min; plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
xline(endpoint - 100);
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize')); legend('ad align', 'targ align')
if save_flag; saveas(gcf, 'dfof align', 'jpg'); end
% close all

range_base = 1:4; range_resp = 14:16;
% prompt = 'base window = 1:3. what is resp window? '; range_resp = input(prompt); close

%% bunnytop early vs late half session resp

ntrial_half = floor(ntrial / 2);

if contains(sess_flag, 'A')
    trial_seq = 1 : ntrial_half;
    dfof_align_ad = dfof_align_ad(:,trial_seq,:);
    id_isi3{1, 3} = id_isi3{1, 3}(id_isi3{1, 3} <= ntrial_half);
    id_ad = id_ad(id_ad <= ntrial_half);
    for i = 1 : nori
        id_ori{i, 1} = id_ori{i, 1}(id_ori{i, 1} <= ntrial_half);
    end
elseif contains(sess_flag, 'B')
    trial_seq = ntrial_half+1 : ntrial;
    dfof_align_ad = dfof_align_ad(:,trial_seq,:);
    id_isi3{1, 3} = id_isi3{1, 3}(id_isi3{1, 3} > ntrial_half);
    id_isi3{1, 3} = id_isi3{1, 3} - ntrial_half;
    id_ad = id_ad(id_ad > ntrial_half);
    id_ad = id_ad - ntrial_half;
    for i = 1 : nori
        id_ori{i, 1} = id_ori{i, 1}(id_ori{i, 1} > ntrial_half);
        id_ori{i, 1} = id_ori{i, 1} - ntrial_half;
    end
end

%% response to adapter & targets. get trace (bunny mode: isi=250 only)
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

%% trial-wise response and baseline

[dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align_ad, save_flag);
dfof_ad_trial = dfof_ad_trial(:,:,3);
dfof_base_trial = dfof_base_trial(:,:,3);

[dfof_tg_trial, dfof_base2_trial] = dfof_resp_trialwise(dfof_align_tg, save_flag);
dfof_tg_trial = dfof_tg_trial(:,:,3);
dfof_base2_trial = dfof_base2_trial(:,:,3);

if save_flag; save resp_base_trialwise.mat dfof_ad_trial dfof_tg_trial dfof_base_trial dfof_base2_trial; end

%% find visually driven cells -> vis_driven.ipynb
