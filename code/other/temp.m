%% get with-adapter 0-deg targ resp

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);
id_delta0 = find(delta_seq == 180);

sig_ttest_tg = pi * ones(ncell, ngap); p_ttest_tg = pi * ones(ncell, ngap);
base_avg_tg = pi * ones(ncell, ngap);
resp_avg_tg = pi * ones(ncell, ngap); resp_ste_tg = pi * ones(ncell, ngap); % standard error 
cp_win_tg = cell(ncell, ngap);
dfof_avg_tg = pi * ones(ncell, ngap); dfof_ste_tg = pi * ones(ncell, ngap); % dF/F

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
       [sig_ttest_tg(icell, igap), p_ttest_tg(icell, igap)] = ttest(base_win{1,igap}, resp_win{1,igap},...
                'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_tg(icell, igap) = mean(base_win{1,igap}); % avg over trials of same ori
        resp_avg_tg(icell, igap) = mean(resp_win{1,igap});
        resp_ste_tg(icell, igap) = std(resp_win{1,igap}) / sqrt(length(resp_win{1,igap}));
        cp_win_tg{icell, igap} = [base_win{1,igap}, resp_win{1,igap}];

        dfof_avg_tg(icell, igap) = mean( (resp_win{1,igap} - base_win{1,igap}) ./ mean(base_win{1,igap}) );
        dfof_ste_tg(icell, igap) = std( (resp_win{1,igap} - base_win{1,igap}) ./ mean(base_win{1,igap}) ) / sqrt(ntrial_cond);
    end
end


%%
sum(sum(sig_ttest_tg,2)>0) % ncells responsive to >= 1 targ ori: 52/103

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

% save with_ad_targ_resp.mat dfof_avg dfof_ste cp_win base_avg resp_avg resp_ste sig_ttest p_ttest
