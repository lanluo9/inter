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

load ori_across_cells.mat
load ori_across_bootstrap_runs_with_replace.mat
load no_ad_targ_resp.mat

%% exp design

ntrial = input.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps
% final trial discarded bc too few frames

cStimOn = double(cell2mat(input.cStimOn)); cStimOff = double(cell2mat(input.cStimOff));
cTarget = celleqel2mat_padded(input.cTargetOn); cTarget = double(cTarget);

adapter_stim_len = min(unique(cStimOff - cStimOn)) % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff) % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input.targetOnTimeMs / frame_rate) % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)) % ITI = 192-193-194 frames, 6.4-6.5 s?
trial_len_list = unique(diff(cStimOn)); % 6.9 or 7.4s

targCon = celleqel2mat_padded(input.tGratingContrast); unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input.tBaseGratingContrast); unique(adapterCon); % adapter contrast 0 or 1

adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
ndelta = length(delta_list); 

[nframe, ncell] = size(npSub_tc) % nframe * ncell

%% align tc by stim onset or targ onset 
% not accounting for ca latency yet. leave room for baseline

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

%% find base & resp window in aligned tc
range = 100; icell = 25;
tt = nanmean(squeeze(tc_trial_align_ad(icell,:,:)), 1); ttt = nanmean(squeeze(tc_trial_align_targ(icell,:,:)), 1);
plot(tt(1:range), 'r'); hold on; plot(ttt(1:range), 'b')
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
legend('ad align', 'targ align')

%% index of trials
ca_latency = 7; % monitor stim onset in frame #1 lead to neural signal in frame #8
isi_list = cTarget - cStimOff;

id_750 = find(isi_list > min(isi_list)); id_750(id_750 > ntrial) = [];
id_250 = find(isi_list == min(isi_list)); id_750(id_750 > ntrial) = [];
id_noad = intersect(find(targCon == 1),find(adapterCon == 0)); id_noad(id_noad > ntrial) = [];
id_ad = intersect(find(targCon == 1),find(adapterCon == 1)); id_ad(id_ad > ntrial) = [];

%% find vis-driven neurons: significant for no-ad targ or ad

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal

sig_ttest = pi * ones(ncell, ndelta); p_ttest = pi * ones(ncell, ndelta);
base_avg = pi * ones(ncell, ndelta);
resp_avg = pi * ones(ncell, ndelta); resp_ste = pi * ones(ncell, ndelta); % standard error 
cp_win = cell(ncell, ndelta);
dfof_avg = pi * ones(ncell, ndelta); dfof_ste = pi * ones(ncell, ndelta); % dF/F

id_noadapter = intersect(find(targCon == 1),find(adapterCon == 0));
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    
    base_win = []; resp_win = [];
    for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
        id_targ = find(targ_start == targ_start_list(igap));
        idx = intersect(intersect(id_targ, id_delta), id_noadapter); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [targ_start_list(igap) - targ_stim_len : targ_start_list(igap) - 1]; % adapted baseline just bef targ onset
        range_targ_resp = [targ_start_list(igap) : targ_start_list(igap) + targ_stim_len - 1]; % targ onset til targ fin

        base_win = [base_win; mean(squeeze(tc_trials(icell, idx, range_adapt_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trials(icell, idx, range_targ_resp)),2)];
    end
        
    [sig_ttest(icell, idelta), p_ttest(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg(icell, idelta) = mean(resp_win);
    resp_ste(icell, idelta) = std(resp_win) / sqrt(length(resp_win));
    cp_win{icell, idelta} = [base_win, resp_win];

    dfof_avg(icell, idelta) = mean( (resp_win - base_win) ./ mean(base_win) );
    dfof_ste(icell, idelta) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrial_cond);
    
end
end

sum(sum(sig_ttest,2)>0) % ncells responsive to >= 1 targ ori: 52/103

subplot(1,2,1)
imagesc(sig_ttest); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest(:,:,1)); colorbar
title('p value')

set(gcf, 'Position', get(0, 'Screensize'));
% cd C:\Users\lan\Documents\repos\inter\code
% saveas(gcf, ['visual_driven_cells_noadapter.jpg'])
% close

% save no_ad_targ_resp.mat dfof_avg dfof_ste cp_win base_avg resp_avg resp_ste sig_ttest p_ttest


