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
load resp_noad_targ.mat
load resp_ad_targ0deg.mat 
load trace_by_cond.mat
load trace_by_cond_dfof.mat

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
% range = 223; 
range = 50;
t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,:,:)), 1));
t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,:,:)), 1)); 
t_tg = squeeze(nanmean(t(:,:), 1)); 

plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b')
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
legend('ad align', 'targ align')

%% index of trials
ca_latency = 7; % monitor stim onset in frame #1 lead to neural signal in frame #8
isi_list = cTarget - cStimOff; 
ngap = length(unique(cTarget - cStimOn));

id_750 = find(isi_list > min(isi_list)); id_750(id_750 > ntrial) = [];
id_250 = find(isi_list == min(isi_list)); id_750(id_750 > ntrial) = [];
id_gaps = {id_250, id_750};
id_noad = intersect(find(targCon == 1),find(adapterCon == 0)); id_noad(id_noad > ntrial) = [];
id_ad = intersect(find(targCon == 1),find(adapterCon == 1)); id_ad(id_ad > ntrial) = [];

%% find vis-driven neurons: significant for no-ad targ or ad

% vis-driven by noad targ
sig_ttest_noad = pi * ones(ncell, ndelta); p_ttest_noad = pi * ones(ncell, ndelta);
base_avg_noad = pi * ones(ncell, ndelta);
resp_avg_noad = pi * ones(ncell, ndelta); resp_ste_noad = pi * ones(ncell, ndelta); % standard error 
cp_win_noad = cell(ncell, ndelta);
dfof_avg_noad = pi * ones(ncell, ndelta); dfof_ste_noad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_noad); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11]; % ref session: find base & resp window in aligned tc

        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2)];
    end
        
    [sig_ttest_noad(icell, idelta), p_ttest_noad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_noad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_noad(icell, idelta) = mean(resp_win);
    resp_ste_noad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));
    cp_win_noad{icell, idelta} = [base_win, resp_win];

    dfof_avg_noad(icell, idelta) = mean( (resp_win - base_win) ./ mean(base_win) );
    dfof_ste_noad(icell, idelta) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrial_cond);
    
end
end
vis_driven_noad = sum(sig_ttest_noad,2)>0;
% save resp_noad_targ.mat dfof_avg_noad dfof_ste_noad cp_win_noad base_avg_noad resp_avg_noad resp_ste_noad sig_ttest_noad p_ttest_noad

sum(sum(sig_ttest_noad,2)>0) % ncells responsive to >= 1 targ ori: 54/103
figure
subplot(1,2,1)
imagesc(sig_ttest_noad); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest_noad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
% saveas(gcf, ['visual_driven_cells_noad_targ.jpg'])
% close


% vis-driven by ad
sig_ttest_ad = pi * ones(ncell, ndelta); p_ttest_ad = pi * ones(ncell, ndelta);
base_avg_ad = pi * ones(ncell, ndelta);
resp_avg_ad = pi * ones(ncell, ndelta); resp_ste_ad = pi * ones(ncell, ndelta); % standard error 
cp_win_ad = cell(ncell, ndelta);
dfof_avg_ad = pi * ones(ncell, ndelta); dfof_ste_ad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11]; % ref session: find base & resp window in aligned tc

        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2)];
    end
        
    [sig_ttest_ad(icell, idelta), p_ttest_ad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_ad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_ad(icell, idelta) = mean(resp_win);
    resp_ste_ad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));
    cp_win_ad{icell, idelta} = [base_win, resp_win];

    dfof_avg_ad(icell, idelta) = mean( (resp_win - base_win) ./ mean(base_win) );
    dfof_ste_ad(icell, idelta) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrial_cond);
    
end
end
vis_driven_ad = sum(sig_ttest_ad,2)>0;
% save resp_ad.mat dfof_avg_ad dfof_ste_ad cp_win_ad base_avg_ad resp_avg_ad resp_ste_ad sig_ttest_ad p_ttest_ad

sum(sum(sig_ttest_ad,2)>0) % ncells responsive to >= 1 targ ori: 69/103
subplot(1,2,1)
imagesc(sig_ttest_ad); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest_ad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
% saveas(gcf, ['visual_driven_cells_adapter.jpg'])
% close

%% cell list by property

vis_driven = vis_driven_ad | vis_driven_noad;
vis_driven_cell_list = find(vis_driven); % vis-driven cells are activated by adapter or no-ad targ

nrun = 1000;
ori_closeness = sort(abs(ori_pref_runs - ori_pref_cells), 2); % sort each row aka cell
percentile_threshold = 0.90;
percentile_idx = percentile_threshold * nrun;
ori_perc = ori_closeness(:, percentile_idx);
good_fit_cell = ori_perc<22.5;
ncell_good_fit = sum(good_fit_cell)

sharp_cell = fit_param(:, 3) > 3;
ncell_sharp = sum(sharp_cell)

best_cell = vis_driven & sharp_cell & good_fit_cell;
best_cell_list = find(best_cell);
ncell_best = length(best_cell_list)

%% resp to with-adapter 0-deg targ 

sig_ttest_tg0 = pi * ones(ncell, ngap); p_ttest_tg0 = pi * ones(ncell, ngap);
base_avg_tg0 = pi * ones(ncell, ngap);
resp_avg_tg0 = pi * ones(ncell, ngap); resp_ste_tg0 = pi * ones(ncell, ngap); % standard error 
cp_win_tg0 = cell(ncell, ngap);
dfof_avg_tg0 = pi * ones(ncell, ngap); dfof_ste_tg0 = pi * ones(ncell, ngap); % dF/F

id_delta0 = find(delta_seq == 180);
for icell = 1 : ncell
    base_win = cell(1,2); resp_win = cell(1,2);
    
    for igap =  1 : ngap
        idx = intersect(intersect(id_gaps{igap}, id_delta0), id_ad); % use only with-adapter 0-deg trials
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11];
        base_win{1,igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_win{1,igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2);
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
% save resp_ad_targ0deg.mat dfof_avg_tg0 dfof_ste_tg0 cp_win_tg0 base_avg_tg0 resp_avg_tg0 resp_ste_tg0 sig_ttest_tg0 p_ttest_tg0

%% Fig 1E: compare ad_resp vs 0-deg targ_resp_250/750 of vis_driven cells
% buggy

pref_0_cell = find(ori_pref_cells > (delta_list(end-1) + delta_list(end))/2 | ori_pref_cells < delta_list(1)/2);
pref_0_cell = pref_0_cell( dfof_avg_tg0(pref_0_cell,1)>0 & dfof_avg_tg0(pref_0_cell,2)>0); % dfof should >0
pref_0_cell = pref_0_cell( dfof_avg_ad(pref_0_cell)>0 ); % dfof should >0
pref_0_cell = intersect(pref_0_cell, vis_driven_cell_list);

[count, idx] = histc(dfof_avg_ad(pref_0_cell),0:0.02:ceil(max(dfof_avg_ad(pref_0_cell))*10)/10);
count_nonzero_id = find(count~=0);
% histogram(dfof_avg_ad(pref_0_cell))
resp_bin_ad = accumarray(idx(:), dfof_avg_ad(pref_0_cell),[],@mean)
resp_bin_tg = [];
resp_bin_std = []; 
resp_bin_ste = zeros(length(count),2);
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@mean);
resp_bin_std(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@std);
resp_bin_std(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@std);
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
% ylim([0,1])
legend('isi 250', 'isi 750')
% saveas(gcf, ['Fig1E cells prefer 0 deg binned by ad-resp.jpg'])
% close

%% Fig 1C (approx): time course of adaptation recovery for all cells
% buggy. smth wrong w dfof_avg_tg0?

isi_sec = [0.250, 0.750, 4];
norm_targ_resp_mean = [mean((dfof_avg_tg0(vis_driven,1)./dfof_avg_ad(vis_driven))), ...
    mean(dfof_avg_tg0(vis_driven,2)./dfof_avg_ad(vis_driven)), 1];

scatter(isi_sec, norm_targ_resp_mean, 'r*'); hold on
f2 = fit(isi_sec', norm_targ_resp_mean','exp1')
plot(f2, 'b')
% xlim([0,4+0.5])
% ylim([0,1+0.3])
xlabel('ISI (s)')
ylabel('normalized dF/F')
legend off
% saveas(gcf, ['Fig1C time course of adaptation recovery of vis-driven cells.jpg'])
% close

%% trial trace by cond or trace noad, converted to dfof

trace_cond = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % with-ad, specific isi & ori
        range_trace = [1 : max(trial_len_list)]; 
        trace_cond{icell, idelta, igap} = squeeze(tc_trial_align_ad(icell, idx, range_trace)); % [ntrial, trial_len]
    end
end
end

trace_no_ad = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_noad); % no-ad, merge isi & specific ori
        range_trace = [1 : max(trial_len_list)]; 
        trace_no_ad{icell, idelta, igap} = squeeze(tc_trial_align_targ(icell, idx, range_trace)); 
    end
end
end


trace_noad_dfof = cell(ncell, ndelta); trace_cond_dfof = cell(ncell, ndelta, ngap);

for icell = 1:ncell
for idelta = 1:ndelta
    trace_base = nanmean(trace_no_ad{icell, idelta}(:, 1:3), 2);
    trace_noad_dfof{icell, idelta} = (trace_no_ad{icell, idelta} - trace_base) ./ trace_base;
end
end

for icell = 1:ncell
for idelta = 1:ndelta
    for igap = 1:ngap
        trace_base = nanmean(trace_cond{icell, idelta, igap}(:, 1:3), 2);
        trace_cond_dfof{icell, idelta, igap} = (trace_cond{icell, idelta, igap} - trace_base) ./ trace_base;
    end
end
end

% cd C:\Users\lan\Documents\repos\inter\code
% save trace_by_cond.mat trace_cond trace_no_ad 
% save trace_by_cond_dfof.mat trace_cond_dfof trace_noad_dfof

%% Fig 2B: cell resp by dir & condition (no_ad, 750, 250) for "best cells"

cell_list_now = vis_driven_cell_list;
for ii = 1 : length(cell_list_now)
icell = cell_list_now(ii);
trace_no_ad_avg = []; trace_cond_avg_750 = []; trace_cond_avg_250 = [];

for idelta = 1 : ndelta
    trace_no_ad_avg(idelta, :) = nanmean(trace_noad_dfof{icell, idelta}, 1);
    trace_cond_avg_750(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 2}, 1);
    trace_cond_avg_250(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 1}, 1);
end

adjusted_trace_len = 3.5 * 30; % trace len 3.5 s according to Jin2019 Fig 2B
trace_no_ad_avg = [trace_no_ad_avg(end,:); trace_no_ad_avg]; % wrap around
trace_no_ad_avg = trace_no_ad_avg(:, 1:adjusted_trace_len); % trim
trace_cond_avg_750 = [trace_cond_avg_750(end,:); trace_cond_avg_750]; 
trace_cond_avg_750 = trace_cond_avg_750(:, 1:adjusted_trace_len);
trace_cond_avg_250 = [trace_cond_avg_250(end,:); trace_cond_avg_250]; 
trace_cond_avg_250 = trace_cond_avg_250(:, 1:adjusted_trace_len);

x = [length(trace_no_ad_avg), length(trace_cond_avg_750), length(trace_cond_avg_250)];
xmax = max(x);
y = [trace_no_ad_avg, trace_cond_avg_750, trace_cond_avg_250];
ymin = min(y(:));
ymax = max(y(:)) + 0.05 * (max(y(:)) - min(y(:)));

figure('units','normalized','outerposition',[0 0 1/2 1]);
for col = 1 : 3
    for idelta = 1 : (ndelta+1)
        subplot(ndelta+1, 3, col + 3*(idelta-1));
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
        yticks(round(ymin*10)/10 : 0.3 : round(ymax*10)/10)
    end
end
saveas(gcf, ['dfof trace ', num2str(icell), '.jpg'])
close    
end