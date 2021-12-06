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

% caiman bunny500 gcamp6s 
dataset_list = struct;
dataset_list.mouse = [1339]; 
dataset_list.date = [210922];
dataset_list.area = {'V1'};

%% load [xls, timecourse, stim]

% for iset = 1 %: length(dataset_list.date)
iset = 1
save_flag = 1; % toggle this to save/skip all .mat creation below

date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset)); imouse = ['i', mouse]
area = dataset_list.area{1,iset}
[~, frame_rate, input_behav, info] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);

areamousedate = [area '_' imouse '_' date];
result_folder = [root_path, '\mat\', areamousedate, '_caiman'];
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

%% substitute npSub_tc w caiman

df = load('C:\Users\ll357\Documents\CaImAn\demos\caiman_activity_i1339_210922_multisess.mat');
df_pile = (df.df)'; % shape = [nframe, ncell]

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
