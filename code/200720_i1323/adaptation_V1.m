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
% saveas(gcf, ['visual_driven_cells_adapter.jpg'])
% close
% save adapter_resp.mat dfof_avg_ad dfof_ste_ad cp_win_ad base_avg_ad resp_avg_ad resp_ste_ad sig_ttest_ad p_ttest_ad

%% get with-adapter 0-deg targ resp

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);
id_delta0 = find(delta_seq == 180);

sig_ttest_tg0 = pi * ones(ncell, ngap); p_ttest_tg0 = pi * ones(ncell, ngap);
base_avg_tg0 = pi * ones(ncell, ngap);
resp_avg_tg0 = pi * ones(ncell, ngap); resp_ste_tg0 = pi * ones(ncell, ngap); % standard error 
cp_win_tg0 = cell(ncell, ngap);
dfof_avg_tg0 = pi * ones(ncell, ngap); dfof_ste_tg0 = pi * ones(ncell, ngap); % dF/F

for icell = 1 : ncell
    base_win = cell(1,2); resp_win = cell(1,2);
    for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
        id_targ = find(targ_start == targ_start_list(igap));
        idx = intersect(intersect(id_targ, id_delta0), id_adapter); % use only with-adapter 0-deg trials
        ntrial_cond = length(idx); 
        
        range_adapt_base = [targ_start_list(igap) - targ_stim_len : targ_start_list(igap) - 1]; % adapted baseline just bef targ onset
        range_targ_resp = [targ_start_list(igap) : targ_start_list(igap) + targ_stim_len - 1]; % targ onset til targ fin
        base_win{1,igap} = mean(squeeze(tc_trials(icell, idx, range_adapt_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_win{1,igap} = mean(squeeze(tc_trials(icell, idx, range_targ_resp)),2);
    end

    for igap =  1 : ngap
       [sig_ttest_tg0(icell, igap), p_ttest_tg0(icell, igap)] = ttest(base_win{1,igap}, resp_win{1,igap},...
                'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_tg0(icell, igap) = mean(base_win{1,igap}); % avg over trials of same ori
        resp_avg_tg0(icell, igap) = mean(resp_win{1,igap});
        resp_ste_tg0(icell, igap) = std(resp_win{1,igap}) / sqrt(length(resp_win{1,igap}));
        cp_win_tg0{icell, igap} = [base_win{1,igap}, resp_win{1,igap}];

        dfof_avg_tg0(icell, igap) = mean( (resp_win{1,igap} - base_win{1,igap}) ./ mean(base_win{1,igap}) );
        dfof_ste_tg0(icell, igap) = std( (resp_win{1,igap} - base_win{1,igap}) ./ mean(base_win{1,igap}) ) / sqrt(ntrial_cond);
    end
end

sum(sig_ttest_tg0, 1) % ncells responsive to with-ad targ whose isi=250 or 750: 0/103 | 5/103

subplot(1,2,1)
imagesc(sig_ttest_tg0); colorbar
title('visually driven by with-adapter targ, isi=250 or 750')
subplot(1,2,2)
imagesc(p_ttest_tg0(:,:,1)); colorbar
title('p value')

set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
% saveas(gcf, ['visual_driven_cells_targ_after_adapter.jpg'])
% close
% save with_ad_0deg_targ_resp.mat dfof_avg_tg0 dfof_ste_tg0 cp_win_tg0 base_avg_tg0 resp_avg_tg0 resp_ste_tg0 sig_ttest_tg0 p_ttest_tg0

%% evaluate best cells

nrun=1000;
ori_closeness = sort(abs(ori_pref_runs - ori_pref_cells), 2); % sort each row aka cell
percentile_threshold = 0.90;
percentile_idx = percentile_threshold * nrun;
ori_perc = ori_closeness(:, percentile_idx);
% sum(ori_perc<22.5) / length(ori_perc)
% histogram(ori_perc,8)

good_fit_cell = ori_perc<22.5;
sum(good_fit_cell)
sig_vis_cell = sum(sig_ttest,2)>0; % sig vis-driven by no-ad targ
sum(sig_vis_cell)
ori_cell = good_fit_cell & sig_vis_cell;
sum(ori_cell)

sharp_cell = fit_param(:, 3) > 3;
real_ori_cell = ori_cell & sharp_cell;
real_ori_cell_list = find(real_ori_cell);
ncell_real_ori = length(real_ori_cell_list)

%% Fig 1E: compare ad_resp vs 0-deg targ_resp_250/750 of vis_driven cells

id_isi_1 = intersect(find(cTarget - cStimOff < 10), id_adapter); % isi 250
id_isi_2 = intersect(find(cTarget - cStimOff >= 10), id_adapter); % isi 750
% use all cells
figure
scatter(dfof_avg_ad, dfof_avg_tg0(:,1)./dfof_avg_ad, 'b*'); hold on % isi 250
scatter(dfof_avg_ad, dfof_avg_tg0(:,2)./dfof_avg_ad, 'r*') % isi 750
yl = ylim; % query [ymin ymax]
% ylim([0, yl(2)])
% ylim([0, 2.5])
xlabel('adapter resp dF/F')
ylabel('0 deg targ resp dF/F norm-ed by ad resp')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['Fig1E all cells.jpg'])
close

% use only 0-preferring cells. other cells' ad & tg resp will certainly co-vary
% pref_0_cell = find(ori_pref_cells > delta_list(end-1));
pref_0_cell = find(ori_pref_cells > (delta_list(end-1) + delta_list(end))/2 | ori_pref_cells < delta_list(1)/2);
pref_0_cell = pref_0_cell( dfof_avg_tg0(pref_0_cell,1)>0 & dfof_avg_tg0(pref_0_cell,2)>0);
% pref_0_cell = find(ori_pref_cells > 170);
figure
subplot(1,2,1)
scatter(dfof_avg_ad(pref_0_cell), dfof_avg_tg0(pref_0_cell,1)./dfof_avg_ad(pref_0_cell), 'b*'); hold on
title('isi=250')
subplot(1,2,2)
scatter(dfof_avg_ad(pref_0_cell), dfof_avg_tg0(pref_0_cell,2)./dfof_avg_ad(pref_0_cell), 'r*')
title('isi=750')
% for n=1:2
% AX_handles(n) = subplot(1,2,n)
% end
% set(AX_handles,'YLim',[-0.04 0.14])
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['Fig1E cells prefer 0 deg.jpg'])
close

% bin all cells by adapter resp
[count_bin, idx] = histc(dfof_avg_ad(vis_driven_ad),0:0.05:ceil(max(dfof_avg_ad(vis_driven_ad))*10)/10);
% histogram(dfof_avg_ad(vis_driven_ad))
resp_bin_ad = accumarray(idx(:),dfof_avg_ad(vis_driven_ad),[],@mean)
resp_bin_tg = [];
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,2),[],@mean)

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b*'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r*')
ylim([0,1])
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['Fig1E all cells binned by ad-resp.jpg'])
close

% bin 0-deg preferring cells by ad resp
[count_bin, idx] = histc(dfof_avg_ad(pref_0_cell), 0:0.05:ceil(max(dfof_avg_ad(pref_0_cell))*10)/10);
% histogram(dfof_avg_ad(pref_0_cell))
resp_bin_ad = accumarray(idx(:),dfof_avg_ad(pref_0_cell),[],@mean)
resp_bin_tg = [];
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@mean)

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b*'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r*')
ylim([0,1])
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['Fig1E cells prefer 0 deg binned by ad-resp.jpg'])
close

%% Fig 1C (approx): time course of adaptation recovery for all cells

% histogram(dfof_avg_tg0(:,1)./dfof_avg_ad)% isi 250
% hold on
% histogram(dfof_avg_tg0(:,2)./dfof_avg_ad)

isi_list = [0.250, 0.750, 4];
norm_targ_resp_list = [median((dfof_avg_tg0(:,1)./dfof_avg_ad)), median(dfof_avg_tg0(:,2)./dfof_avg_ad), 1];
scatter(isi_list, norm_targ_resp_list)
hold on
f = fit(isi_list',norm_targ_resp_list','exp1')
plot(f,isi_list',norm_targ_resp_list')
xlim([0,4])
ylim([0,1])
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['Fig1C time course of adaptation recovery w only 2 data points.jpg'])
close


%% Fig 2B: cell resp by dir & condition (no_ad, 750, 250)

for icell = 1:nell
    
    
    
end
