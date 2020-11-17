function vis_cell_criteria(sig_alpha)
global nori ori_seq ori_list ncell nisi id_isi

% vis-driven by any noad targ
sig_ttest = pi * ones(ncell, ndelta); p_ttest = pi * ones(ncell, ndelta);
base_avg = pi * ones(ncell, ndelta); resp_avg = pi * ones(ncell, ndelta); 
dfof_avg = pi * ones(ncell, ndelta); dfof_ste = pi * ones(ncell, ndelta); % dF/F

for iori = 1 : nori
    id_ori = find(ori_seq == ori_list(iori));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for iisi =  1 : nisi % order: 750, 250
        idx = intersect(intersect(id_isi{iisi}, id_ori), id_noad); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_resp)),2)];
    end
        
    [sig_ttest(icell, iori), p_ttest(icell, iori)] = ttest(base_win, resp_win,...
            'alpha',sig_alpha./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg(icell, iori) = mean(base_win); % avg over trials of same ori
    resp_avg(icell, iori) = mean(resp_win);
    resp_ste(icell, iori) = std(resp_win) / sqrt(length(resp_win));

    dfof_avg(icell, iori) = mean( resp_win - base_win );
    dfof_ste(icell, iori) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    
end
end
