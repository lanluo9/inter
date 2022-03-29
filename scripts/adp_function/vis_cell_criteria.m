function [sig_vis, p_vis, dfof_avg] = vis_cell_criteria(dfof_align, compared_obj, sig_alpha)

% cells visually driven by
% adapter (0 deg): one-way t-test response vs baseline
% or any noad target (8 oris): one-way t-test w bonferonni correction for number of stimuli
% 
% input: dfof_align by ad or tg, linked w/ compared_obj = 'ad' or 'tg_any'. sig_alpha = sig threshold
% output: sig_vis & p_vis. dfof_avg for ad or noad-tg (all oris)

global nori id_ori id_noad id_ad range_base range_resp ncell

switch compared_obj
case 'ad' % vis-driven by ad: 1-way t-test
    for icell = 1 : ncell
        base_win = []; resp_win = [];
        base_win = [base_win; mean(squeeze(dfof_align(icell, id_ad, range_base)),2)]; 
        resp_win = [resp_win; mean(squeeze(dfof_align(icell, id_ad, range_resp)),2)];

        [sig_vis(icell), p_vis(icell)] = ttest(base_win, resp_win,...
                'alpha',sig_alpha, 'tail','left'); 
        base_avg(icell) = mean(base_win); resp_avg(icell) = mean(resp_win);
        dfof_avg(icell) = mean( resp_win - base_win );
    end

case 'tg_any' % vis-driven by any noad targ: 1-way t-test w bonferonni correction
    for icell = 1 : ncell
        base_win = []; resp_win = [];

    for iori = 1 : nori
        idx = intersect(id_ori{iori}, id_noad); % noad trials with 1 ISI & 1 ori
%         ntrial_cond = length(idx); 
        base_win = [base_win; mean(squeeze(dfof_align(icell, idx, range_base)),2)]; 
        resp_win = [resp_win; mean(squeeze(dfof_align(icell, idx, range_resp)),2)];

        [sig_vis(icell, iori), p_vis(icell, iori)] = ttest(base_win, resp_win,...
                'alpha',sig_alpha./(nori - 1), 'tail','left'); % sig = base<resp, Bonferroni correction for nori
        base_avg(icell, iori) = mean(base_win); resp_avg(icell, iori) = mean(resp_win);
        dfof_avg(icell, iori) = mean( resp_win - base_win );
    %     dfof_ste(icell, iori) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    end
    end

end
