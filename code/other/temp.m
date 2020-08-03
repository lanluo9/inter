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

%%

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

%%
sum(sig_ttest_cond, 1) % ncells responsive to with-ad targ whose isi=250 or 750: 0/103 | 5/103

subplot(1,2,1)
imagesc(sig_ttest_cond); colorbar
title('visually driven by with-adapter targ, isi=250 or 750')
subplot(1,2,2)
imagesc(p_ttest_cond(:,:,1)); colorbar
title('p value')

set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
% saveas(gcf, ['visual_driven_cells_targ_with_adapter.jpg'])
% close
save with_ad_all_oris_targ_resp.mat dfof_avg_cond dfof_ste_cond cp_win_cond base_avg_cond resp_avg_cond resp_ste_cond sig_ttest_cond p_ttest_cond
