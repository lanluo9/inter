function [sig_vis, p_vis] = vis_driven_ohki(dfof_align, sig_alpha)

% cells visually driven by any noad target: 1 way ANOVA across [8 ori + 1 baseline]

global nori id_ori id_noad range_base range_resp ncell

for icell = 1 : ncell
    base_win = []; resp_win = [];

for iori = 1 : nori
    idx = intersect(id_ori{iori}, id_noad); % noad trials with 1 ISI & 1 ori
    base_win = [base_win; mean(squeeze(dfof_align(icell, idx, range_base)),2)]; 
    resp_win = [resp_win; mean(squeeze(dfof_align(icell, idx, range_resp)),2)];

%     [sig_vis(icell, iori), p_vis(icell, iori)] = ttest(base_win, resp_win,...
%             'alpha',sig_alpha./(nori - 1), 'tail','left'); % sig = base<resp, Bonferroni correction for nori

%     base_avg(icell, iori) = mean(base_win); resp_avg(icell, iori) = mean(resp_win);
%     dfof_avg(icell, iori) = mean( resp_win - base_win );
%     dfof_ste(icell, iori) = std( resp_win - base_win ) / sqrt(ntrial_cond);
end
end

end
