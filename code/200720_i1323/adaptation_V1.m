%% set up directory

close all
clear
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

% need automation here. import excel
date = '200720';
ImgFolder = '003';
time = strvcat('1133'); 
mouse = 'i1323';
frame_rate = 30;

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
fnout = fullfile(ll_fn, 'Analysis\2P\', datemouserun); 

%% load data

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"

CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

mouse = '1323';
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp

cd C:\Users\lan\Documents\repos\inter\code\200720_i1323

%% exp design

ntrial = input.trialSinceReset; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps

cStimOn = cell2mat(input.cStimOn);
cStimOff = cell2mat(input.cStimOff);
cTarget = celleqel2mat_padded(input.cTargetOn); cTarget = int64(cTarget);

adapter_stim_len = min(unique(cStimOff - cStimOn)) % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff) % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input.targetOnTimeMs / frame_rate) % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)) % ITI = 192-193-194 frames, 6.4-6.5 s?
trial_len_list = unique(diff(cStimOn)); % 6.9 or 7.4s

targCon = celleqel2mat_padded(input.tGratingContrast);
unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input.tBaseGratingContrast);
unique(adapterCon); % adapter contrast 0 or 1

adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
ndelta = length(delta_list); 

[nframe, ncell] = size(npSub_tc) % nframe * ncell

%% align by trial start

trial_len = diff(cStimOn);
trial_len_final = nframe - cStimOn(end);
trial_len = [trial_len, trial_len_final];

tc_trials = zeros(ncell, ntrial, max(trial_len));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial);
        tc_trials(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len(itrial) - 1); ...
            NaN(max(trial_len) - trial_len(itrial), 1)];
    end
end

%% load target response
ca_latency = 7; % Ca signal latency: signal from frame #1 shows up at frame #8
ad_start = 1 + ca_latency;
id_noadapter = intersect(find(targCon == 1),find(adapterCon == 0));
id_adapter = intersect(find(targCon == 1),find(adapterCon == 1));

% load no_ad_targ_resp.mat

%% get adapter response
clc

sig_ttest_ad = pi * ones(ncell,1); p_ttest_ad = pi * ones(ncell,1);
base_avg_ad = pi * ones(ncell,1);
resp_avg_ad = pi * ones(ncell,1); resp_ste_ad = pi * ones(ncell,1); % standard error 
cp_win_ad = cell(ncell,1);
dfof_avg_ad = pi * ones(ncell,1); dfof_ste_ad = pi * ones(ncell,1); % dF/F
    
for icell = 1 : ncell
    
    base_win = []; resp_win = [];
    idx = id_adapter; % use only with-adapter trials with 1 ISI & 1 ori
    ntrial_cond = length(idx); 

    range_base = [1 : ad_start - 1]; % baseline just bef ad onset
    range_ad_resp = [ad_start : ad_start + adapter_stim_len - 1]; % ad onset til ad fin
    base_win = [base_win; mean(squeeze(tc_trials(icell, idx, range_base)),2)]; % avg over window -> [ntrial_ori, 1]
    resp_win = [resp_win; mean(squeeze(tc_trials(icell, idx, range_ad_resp)),2)];
        
    [sig_ttest_ad(icell), p_ttest_ad(icell)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_ad(icell) = mean(base_win); % avg over trials of same ori
    resp_avg_ad(icell) = mean(resp_win);
    resp_ste_ad(icell) = std(resp_win) / sqrt(length(resp_win));
    cp_win_ad{icell} = [base_win, resp_win];

    dfof_avg_ad(icell) = mean( (resp_win - base_win) ./ mean(base_win) );
    dfof_ste_ad(icell) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrial_cond);
    
end

sum(sum(sig_ttest_ad,2)>0) % ncells responsive to adapter: 77/103 (52/103 estimated by no-adapter targ resp)
vis_driven_ad = sum(sig_ttest_ad,2)>0;
histogram(dfof_avg_ad(vis_driven_ad))

subplot(1,2,1)
imagesc(sig_ttest_ad); colorbar
title('visually driven by adapter')
subplot(1,2,2)
imagesc(p_ttest_ad(:,:,1)); colorbar
title('p value')

set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
saveas(gcf, ['visual_driven_cells_adapter.jpg'])
% close
% save adapter_resp.mat dfof_avg_ad dfof_ste_ad cp_win_ad base_avg_ad resp_avg_ad resp_ste_ad sig_ttest_ad p_ttest_ad

%% compare ad_resp vs 0-deg targ_resp_250/750 of vis_driven cells

id_isi_1 = intersect(find(cTarget - cStimOff < 10), id_adapter); % isi 250
id_isi_2 = intersect(find(cTarget - cStimOff >= 10), id_adapter); % isi 750

% dfof_avg_ad
% bin cells by adapter resp