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

load adapter_resp.mat
vis_driven_ad = sum(sig_ttest_ad,2)>0;

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

load with_ad_0deg_targ_resp.mat

%% evaluate best cells

nrun = 1000;
ori_closeness = sort(abs(ori_pref_runs - ori_pref_cells), 2); % sort each row aka cell
percentile_threshold = 0.90;
percentile_idx = percentile_threshold * nrun;
ori_perc = ori_closeness(:, percentile_idx);
% sum(ori_perc<22.5) / length(ori_perc)
% histogram(ori_perc,8)

good_fit_cell = ori_perc<22.5;
ncell_good_fit = sum(good_fit_cell)
sig_vis_cell = sum(sig_ttest,2)>0; % sig vis-driven by no-ad targ
ncell_sig_vis = sum(sig_vis_cell)
ori_cell = good_fit_cell & sig_vis_cell;
ncell_ori = sum(ori_cell)

sharp_cell = fit_param(:, 3) > 3;
best_cell = ori_cell & sharp_cell;
best_cell_list = find(best_cell);
ncell_best = length(best_cell_list)

%% Fig 1E: compare ad_resp vs 0-deg targ_resp_250/750 of vis_driven cells

% id_isi_1 = intersect(find(cTarget - cStimOff < 10), id_adapter); % isi 250
% id_isi_2 = intersect(find(cTarget - cStimOff >= 10), id_adapter); % isi 750

% %%% use all cells
% figure
% scatter(dfof_avg_ad, dfof_avg_tg0(:,1)./dfof_avg_ad, 'b*'); hold on % isi 250
% scatter(dfof_avg_ad, dfof_avg_tg0(:,2)./dfof_avg_ad, 'r*') % isi 750
% yl = ylim; % query [ymin ymax]
% % ylim([0, yl(2)])
% % ylim([0, 2.5])
% xlabel('adapter resp dF/F')
% ylabel('0 deg targ resp dF/F norm-ed by ad resp')
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1E all cells.jpg'])
% close

%%% use only 0-preferring cells. other cells' ad & tg resp will certainly co-vary
% pref_0_cell = find(ori_pref_cells > delta_list(end-1));
pref_0_cell = find(ori_pref_cells > (delta_list(end-1) + delta_list(end))/2 | ori_pref_cells < delta_list(1)/2);
pref_0_cell = pref_0_cell( dfof_avg_tg0(pref_0_cell,1)>0 & dfof_avg_tg0(pref_0_cell,2)>0); % dfof should >0
% pref_0_cell = find(ori_pref_cells > 170);

% figure
% subplot(1,2,1)
% scatter(dfof_avg_ad(pref_0_cell), dfof_avg_tg0(pref_0_cell,1)./dfof_avg_ad(pref_0_cell), 'b*'); hold on
% title('isi=250')
% subplot(1,2,2)
% scatter(dfof_avg_ad(pref_0_cell), dfof_avg_tg0(pref_0_cell,2)./dfof_avg_ad(pref_0_cell), 'r*')
% title('isi=750')
% % for n=1:2
% % AX_handles(n) = subplot(1,2,n)
% % end
% % set(AX_handles,'YLim',[-0.04 0.14])
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1E cells prefer 0 deg.jpg'])
% close

%%% bin all cells by adapter resp
[count, idx] = histc(dfof_avg_ad(vis_driven_ad),0:0.05:ceil(max(dfof_avg_ad(vis_driven_ad))*10)/10);
count = count(1:end-1);
count_nonzero_id = find(count~=0);
% histogram(dfof_avg_ad(vis_driven_ad))
resp_bin_ad = accumarray(idx(:),dfof_avg_ad(vis_driven_ad),[],@mean)
resp_bin_tg = [];
resp_bin_std = []; 
resp_bin_ste = zeros(length(count),2);
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,2),[],@mean)
resp_bin_std(:,1) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,1),[],@std);
resp_bin_std(:,2) = accumarray(idx(:),dfof_avg_tg0(vis_driven_ad,2),[],@std)
resp_bin_ste(count~=0, :) = resp_bin_std(count~=0, :) ./ sqrt(count(count~=0)) % use std or ste?

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b.'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r.')
errorbar(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, resp_bin_std(:,1), 'b', 'LineStyle','none');
errorbar(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, resp_bin_std(:,2), 'r', 'LineStyle','none');
for itext = 1 : length(count_nonzero_id)
    text(resp_bin_ad(count_nonzero_id(itext)), ...
        resp_bin_tg(count_nonzero_id(itext),2)./resp_bin_ad(count_nonzero_id(itext)) + resp_bin_std(count_nonzero_id(itext),2) + 0.02, ...
        ['n=', num2str(count(count_nonzero_id(itext)))], 'HorizontalAlignment', 'center')
end
ylim([0,1])
legend('isi 250', 'isi 750')
set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1E all cells binned by ad-resp.jpg'])
% close

%%% bin 0-deg preferring cells by ad resp
[count, idx] = histc(dfof_avg_ad(pref_0_cell),0:0.05:ceil(max(dfof_avg_ad(pref_0_cell))*10)/10);
count = count(1:end-1);
count_nonzero_id = find(count~=0);
% histogram(dfof_avg_ad(pref_0_cell))
resp_bin_ad = accumarray(idx(:),dfof_avg_ad(pref_0_cell),[],@mean)
resp_bin_tg = [];
resp_bin_std = []; 
resp_bin_ste = zeros(length(count),2);
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@mean)
resp_bin_std(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@std);
resp_bin_std(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@std)
resp_bin_ste(count~=0, :) = resp_bin_std(count~=0, :) ./ sqrt(count(count~=0)) % use std or ste?

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b.'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r.')
errorbar(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, resp_bin_std(:,1), 'b', 'LineStyle','none');
errorbar(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, resp_bin_std(:,2), 'r', 'LineStyle','none');
for itext = 1 : length(count_nonzero_id)
    text(resp_bin_ad(count_nonzero_id(itext)), ...
        resp_bin_tg(count_nonzero_id(itext),2)./resp_bin_ad(count_nonzero_id(itext)) + resp_bin_std(count_nonzero_id(itext),2) + 0.02, ...
        ['n=', num2str(count(count_nonzero_id(itext)))], 'HorizontalAlignment', 'center')
end
ylim([0,1])
legend('isi 250', 'isi 750')
set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1E cells prefer 0 deg binned by ad-resp.jpg'])
% close

%% Fig 1C (approx): time course of adaptation recovery for all cells

vis_driven_ad_pos = dfof_avg_tg0(:,1)>0 & dfof_avg_tg0(:,2)>0 & vis_driven_ad;
sum(vis_driven_ad_pos)
% histogram(dfof_avg_tg0(vis_driven_ad_pos,1)./dfof_avg_ad(vis_driven_ad_pos), 10)% isi 250
% hold on
% histogram(dfof_avg_tg0(vis_driven_ad_pos,2)./dfof_avg_ad(vis_driven_ad_pos), 10)

isi_list = [0.250, 0.750, 4];
% norm_targ_resp_median = [median((dfof_avg_tg0(vis_driven_ad,1)./dfof_avg_ad(vis_driven_ad))), ...
%     median(dfof_avg_tg0(vis_driven_ad,2)./dfof_avg_ad(vis_driven_ad)), ...
%     1];
norm_targ_resp_mean = [mean((dfof_avg_tg0(vis_driven_ad,1)./dfof_avg_ad(vis_driven_ad))), ...
    mean(dfof_avg_tg0(vis_driven_ad,2)./dfof_avg_ad(vis_driven_ad)), ...
    1];

hold on
% scatter(isi_list, norm_targ_resp_median, 'b*')
% f1 = fit(isi_list', norm_targ_resp_median','exp1')
% plot(f1, 'b')
scatter(isi_list, norm_targ_resp_mean, 'r*')
f2 = fit(isi_list', norm_targ_resp_mean','exp1')
plot(f2, 'r')
xlim([0,4+0.5])
ylim([0,1+0.1])
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1C time course of adapstation recovery of vis-driven cells.jpg'])
% close

%% get with-adapter targ resp by dir & isi

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);

base_cond = cell(ncell, ndelta, ngap); resp_cond = cell(ncell, ndelta, ngap);
for icell = 1 : ncell
    
for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
        id_targ = find(targ_start == targ_start_list(igap));
        idx = intersect(intersect(id_targ, id_delta), id_adapter); % with-ad, specific isi & ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [targ_start_list(igap) - targ_stim_len : targ_start_list(igap) - 1]; % adapted baseline just bef targ onset
        range_targ_resp = [targ_start_list(igap) : targ_start_list(igap) + targ_stim_len - 1]; % targ onset til targ fin
        base_cond{icell, idelta, igap} = mean(squeeze(tc_trials(icell, idx, range_adapt_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_cond{icell, idelta, igap} = mean(squeeze(tc_trials(icell, idx, range_targ_resp)),2);
    end
end
end

sig_ttest_cond = pi * ones(ncell, ndelta, ngap); p_ttest_cond = pi * ones(ncell, ndelta, ngap);
base_avg_cond = pi * ones(ncell, ndelta, ngap);
resp_avg_cond = pi * ones(ncell, ndelta, ngap); resp_ste_cond = pi * ones(ncell, ndelta, ngap); % standard error 
cp_win_cond = cell(ncell, ndelta, ngap);
dfof_avg_cond = pi * ones(ncell, ndelta, ngap); dfof_ste_cond = pi * ones(ncell, ndelta, ngap); % dF/F

for icell = 1 : ncell
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap
       [sig_ttest_cond(icell, idelta, igap), p_ttest_cond(icell, idelta, igap)] = ttest(base_cond{icell, idelta, igap}, resp_cond{icell, idelta, igap},...
                'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_cond(icell, idelta, igap) = mean(base_cond{icell, idelta, igap}); % avg over trials of same ori
        resp_avg_cond(icell, idelta, igap) = mean(resp_cond{icell, idelta, igap});
        resp_ste_cond(icell, idelta, igap) = std(resp_cond{icell, idelta, igap}) / sqrt(length(resp_cond{icell, idelta, igap}));
        cp_win_cond{icell, idelta, igap} = [base_cond{icell, idelta, igap}, resp_cond{icell, idelta, igap}];

        dfof_avg_cond(icell, idelta, igap) = mean( (resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap}) ./ mean(base_cond{icell, idelta, igap}) );
        dfof_ste_cond(icell, idelta, igap) = std( (resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap}) ./ mean(base_cond{icell, idelta, igap}) ) / sqrt(ntrial_cond);
    end
end
end

% cd C:\Users\lan\Documents\repos\inter\code
% save with_ad_all_oris_targ_resp.mat dfof_avg_cond dfof_ste_cond cp_win_cond base_avg_cond resp_avg_cond resp_ste_cond sig_ttest_cond p_ttest_cond

load with_ad_all_oris_targ_resp.mat

%% get all cell trial trace by dir & isi

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);

trace_cond = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        id_targ = find(targ_start == targ_start_list(igap));
        
        idx = intersect(intersect(id_targ, id_delta), id_adapter); % with-ad, specific isi & ori
        range_trace = [1 : max(trial_len)]; 
        trace_cond{icell, idelta, igap} = squeeze(tc_trials(icell, idx, range_trace)); % [ntrial, trial_len]
    end
end
end

trace_no_ad = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        id_targ = find(targ_start == targ_start_list(igap));
        idx = intersect(intersect(id_targ, id_delta), id_noadapter); % no-ad, merge isi & specific ori
        range_trace = [1 : max(trial_len)]; 
        trace_no_ad{icell, idelta, igap} = squeeze(tc_trials(icell, idx, range_trace)); 
    end
end
end

cd C:\Users\lan\Documents\repos\inter\code
% save trace_by_cond.mat trace_cond trace_no_ad

load trace_by_cond.mat

%% align no-ad targ resp across fake isi

alt_targ_start = max(targ_start_list) - (min(targ_start_list) - 1);
trace_no_ad_align = cell(ncell, ndelta, ngap); 

for icell = 1 : ncell    
for idelta = 1 : ndelta 
    
    for igap =  1 : ngap 
        if igap == 1
            range_trace = [1 : (1 + max(trial_len) - alt_targ_start)]; % take 1:208 of isi_250 trials
            trace_no_ad_align{icell, idelta, igap} = trace_no_ad{icell, idelta, igap}(:,range_trace);
        elseif igap == 2
            range_trace = [alt_targ_start : max(trial_len)]; % take 16:end of isi_750 trials
            trace_no_ad_align{icell, idelta, igap} = trace_no_ad{icell, idelta, igap}(:,range_trace);
        end
    end
end
end

trace_no_ad_merge = cell(ncell, ndelta);
for icell = 1 : ncell    
for idelta = 1 : ndelta
    trace_no_ad_merge{icell, idelta} = [trace_no_ad_align{icell, idelta, 1}; trace_no_ad_align{icell, idelta, 2}];
end
end

% save('trace_by_cond.mat', 'trace_no_ad_merge', '-append') 

%% convert trace to dfof value

base_range = [80:200]; % according to find_ca_latency.jpg
trace_no_ad_merge_dfof = cell(ncell, ndelta);
trace_cond_dfof = cell(ncell, ndelta, ngap);

for icell = 1:ncell
for idelta = 1:ndelta
    trace_base = nanmean(trace_no_ad_merge{icell, idelta}(:,base_range), 2);
    trace_no_ad_merge_dfof{icell, idelta} = (trace_no_ad_merge{icell, idelta} - trace_base) ./ trace_base;
end
end

for icell = 1:ncell
for idelta = 1:ndelta
    for igap = 1:ngap
        trace_base = nanmean(trace_cond{icell, idelta, igap}(:,base_range), 2);
        trace_cond_dfof{icell, idelta, igap} = (trace_cond{icell, idelta, igap} - trace_base) ./ trace_base;
    end
end
end

%% Fig 2B: cell resp by dir & condition (no_ad, 750, 250) for "best cells"

for ibestcell = 1:length(best_cell_list)
icell = best_cell_list(ibestcell);
trace_no_ad_avg = []; trace_cond_avg_750 = []; trace_cond_avg_250 = [];

for idelta = 1 : ndelta
    trace_no_ad_avg(idelta, :) = nanmean(trace_no_ad_merge_dfof{icell, idelta}, 1);
    trace_cond_avg_750(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 2}, 1);
    trace_cond_avg_250(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 1}, 1);
end
x = [length(trace_no_ad_avg), length(trace_cond_avg_750), length(trace_cond_avg_250)];
xmax = max(x);
y = [trace_no_ad_avg, trace_cond_avg_750, trace_cond_avg_250];
ymin = min(y(:));
ymax = max(y(:));

figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle_LL(num2str(icell))
for col = 1 : 3
    for idelta = 1 : ndelta
        subplot(ndelta, 3, col + 3*(idelta-1));
        if col == 1, plot(trace_no_ad_avg(idelta, :))
            if idelta == 1, title('no adapter'), end
        elseif col == 2, plot(trace_cond_avg_750(idelta, :))
             if idelta == 1, title('isi 750'), end
        elseif col == 3, plot(trace_cond_avg_250(idelta, :))
             if idelta == 1, title('isi 250'), end
        end
%         axis off
        xlim([0, xmax])
        ylim([ymin, ymax])
        xticks(0 : 30 : xmax)
        yticks(round(ymin*10)/10 : 0.1 : round(ymax*10)/10)
    end
end
saveas(gcf, ['dfof trace ', num2str(icell), '.jpg'])
close
    
end

%%
