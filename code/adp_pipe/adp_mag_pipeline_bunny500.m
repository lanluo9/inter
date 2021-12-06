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
% [~, frame_rate, input_behav, info] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);

% get stim mat: input_behav for each session
xls_dir = fullfile(data_fn, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx');
data_now_meta = readtable(xls_file.name);
frame_rate = data_now_meta.(5)(end);
bunny500_id = find(contains(data_now_meta.StimulusSet,'bunny_randStim1_doSameStim2'))

clear input_behav_seq
for i = 1:length(bunny500_id)
    id = bunny500_id(i);
    time = data_now_meta.(8)(id);
    ImgFolder = data_now_meta.(1){id}(1:3);

    fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' num2str(time) '.mat']);
    temp = load(fName); % load behavior data "input"
    input_behav_seq(i) = temp.input; clear temp
end

areamousedate = [area '_' imouse '_' date];
result_folder = [root_path, '\mat\', areamousedate, '_caiman'];
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder)

%% substitute npSub_tc w caiman

df = load('C:\Users\ll357\Documents\CaImAn\demos\caiman_activity_i1339_210922_multisess.mat');
df_pile = (df.df)'; % need to reshape to [nframe_sum, ncell]

t = cellfun(@size,df_pile,'uni',false); 
ncell = size(t,2);
t = cell2mat(t(:,1));
nframe_seq = t(:,2);

for icell = 1:ncell
    
end

%% concat trial stim info for each session
% index by adapter contrast, target ori, isi

ntrial = 0;
adapter_id = [];
target_id = [];
frame_ad = [];
frame_ad_off = [];
frame_tg = [];

for i = 1 : length(bunny500_id)
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
nrep_stim = unique(t(:,2)) % bunny500 3sess gives 3/4/5 rep of each img
