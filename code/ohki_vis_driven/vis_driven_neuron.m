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

% for iset = 1 : length(dataset_list.date)
iset = 1;
date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset))
area = dataset_list.area{1,iset}

imouse = ['i', mouse];
xls_dir = fullfile(data_fn, imouse, date);
cd(xls_dir)
xls_file = dir('*.xlsx'); 
clear dataset_meta
dataset_meta = readtable(xls_file.name); 
time = dataset_meta.(8){end}
ImgFolder = dataset_meta.(1){end}(1:3);
frame_rate = str2num(dataset_meta.(5){end});

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' imouse];
datemouserun = [date '_' imouse '_' run_str];
areamousedate = [dataset_list.area{1,iset} '_' imouse '_' date];

fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"
input_behav = input;
clear input

CD = fullfile(data_fn, imouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp
% fix bug later: retinotopy folder in Analysis naming convention should adhere to tc folder

result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
result_folder = fullfile(result_prefix, areamousedate);
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder);

%% dataset params

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps % final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc); % nframe * ncell

cStimOn = double(cell2mat(input_behav.cStimOn)); cStimOff = double(cell2mat(input_behav.cStimOff));
cTarget = celleqel2mat_padded(input_behav.cTargetOn); cTarget = double(cTarget);

adapter_stim_len = min(unique(cStimOff - cStimOn)); % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff); % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input_behav.targetOnTimeMs / frame_rate); % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)); % ITI = 192-193-194 frames, 6.4-6.5 s?
trial_len_list = unique(diff(cStimOn)); % 6.9 or 7.4s

targCon = celleqel2mat_padded(input_behav.tGratingContrast); unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input_behav.tBaseGratingContrast); unique(adapterCon); % adapter contrast 0 or 1

adapterDir = celleqel2mat_padded(input_behav.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
% fix bug later: need to convert 180 to 0
ndelta = length(delta_list); 

%% index of trials

isi_list = cTarget - cStimOff; 
ngap = length(unique(cTarget - cStimOn));

id_750 = find(isi_list > mean(isi_list)); id_750(id_750 > ntrial) = [];
id_250 = find(isi_list < mean(isi_list)); id_250(id_250 > ntrial) = [];
id_gaps = {id_750, id_250}; 
id_noad = intersect(find(targCon == 1),find(adapterCon == 0)); id_noad(id_noad > ntrial) = [];
id_ad = intersect(find(targCon == 1),find(adapterCon == 1)); id_ad(id_ad > ntrial) = [];

%% align tc by adapter onset or targ onset. normalize aligned tc by 1-sec-long "trial baseline"

trial_len_ad = diff(cStimOn);
tc_trial_align_ad = zeros(ncell, ntrial, max(trial_len_ad));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial);
        tc_trial_align_ad(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len_ad(itrial) - 1); ...
            NaN(max(trial_len_ad) - trial_len_ad(itrial), 1)];
    end
end

trial_len_targ = diff(cTarget);
tc_trial_align_targ = zeros(ncell, ntrial, max(trial_len_targ));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cTarget(itrial);
        tc_trial_align_targ(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len_targ(itrial) - 1); ...
            NaN(max(trial_len_targ) - trial_len_targ(itrial), 1)];
    end
end

trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
tc_trial_base = zeros(ncell, ntrial, trial_base_len);
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial) - trial_base_len;
        tc_trial_base(icell, itrial, :) = [npSub_tc_cell(cStimOn(itrial) - trial_base_len : cStimOn(itrial) - 1)];
    end
end
t = squeeze(nanmean(squeeze(tc_trial_base(:,:,:)), 1)); 
t_base = squeeze(nanmean(t(:,:), 1)); 
plot(t_base, 'k');

result_sub = fullfile(result_folder, 'pre-processing');
if ~exist(result_sub); mkdir(result_sub); end
cd(result_sub)
% saveas(gcf, ['trial base'], 'jpg'); close 

tc_trial_base_avg = mean(tc_trial_base, 3);
for icell = 1:ncell
    for itrial = 1:ntrial
        tc_trial_align_ad(icell, itrial, :) = tc_trial_align_ad(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) -1;
        tc_trial_align_targ(icell, itrial, :) = tc_trial_align_targ(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) -1;
    end
end

%% find base & resp window in aligned tc
% range = 223; 
range = 50;
t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,:,:)), 1));
t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,:,:)), 1)); 
t_tg = squeeze(nanmean(t(:,:), 1)); 

plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
legend('ad align', 'targ align')
saveas(gcf, ['aligned_tc_zoomin'], 'jpg')
close

% t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,:,:)), 1));
% t_ad = squeeze(nanmean(t(:,:), 1)); 
% t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_750,:)), 1)); 
% t_tg_750 = squeeze(nanmean(t(:,:), 1)); 
% t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_250,:)), 1)); 
% t_tg_250 = squeeze(nanmean(t(:,:), 1)); 
% 
% plot(t_ad(1:range), 'r'); hold on; plot(t_tg_750(1:range), 'b'); plot(t_tg_250(1:range), 'g'); 
% grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
% legend('ad align', 'targ align 750', 'targ align 250')

% % t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_750,:)), 1));
% % t_ad_750 = squeeze(nanmean(t(:,:), 1)); 
% % t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_250,:)), 1));
% % t_ad_250 = squeeze(nanmean(t(:,:), 1)); 
% % 
% % t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_750,:)), 1)); 
% % t_tg_750 = squeeze(nanmean(t(:,:), 1)); 
% % t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_250,:)), 1)); 
% % t_tg_250 = squeeze(nanmean(t(:,:), 1)); 
% % 
% % plot(t_ad_750(1:range), 'r'); hold on; plot(t_ad_250(1:range), 'r--'); 
% % plot(t_tg_750(1:range), 'b'); plot(t_tg_250(1:range), 'g'); 
% % grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
% % legend('ad align 750', 'ad align 250', 'targ align 750', 'targ align 250')


range_base = [1:3]; 
range_resp = [9:12];
% prompt = 'base window = 1:3. what is resp window? ';
% range_resp = input(prompt)
% close

%% find vis-driven & good-fit