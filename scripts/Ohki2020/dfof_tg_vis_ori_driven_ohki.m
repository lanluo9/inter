function [dfof_tg_avg, dfof_tg_sem, dfof_tg_std, sig_ori_driven, sig_vis_driven] = dfof_tg_vis_ori_driven_ohki(dfof_align, sig_alpha, resp_threshold)

% output: dfof avg & sem & std for tg = ncell x nori x nisi3 [noad 750 250]

global ncell nori nisi id_isi3 id_ori range_base range_resp

base_cond = cell(ncell, nori, nisi); resp_cond = cell(ncell, nori, nisi);
init = pi * ones(ncell, nori, nisi); dfof_tg_avg = init; dfof_tg_sem = init; dfof_tg_std = init;
sig_ori_driven = pi * ones(ncell, nori);
sig_vis_driven = zeros(ncell,1);

for icell = 1 : ncell
    baseline = []; ori_resp = []; ori_list = [];
    
    for iori = 1 : nori 
        for iisi =  1 : length(id_isi3) 
            idx = intersect(id_ori{iori}, id_isi3{iisi}); 
            base_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_base)),2); 
            resp_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_resp)),2);

            ntrial_cond = length(base_cond{icell, iori, iisi});
            dfof_tg_avg(icell, iori, iisi) = mean( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
            dfof_tg_sem(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} ) / sqrt(ntrial_cond);
            dfof_tg_std(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
        end
        
        [sig_ori_driven(icell, iori), ~] = ttest(base_cond{icell, iori, 1}, resp_cond{icell, iori, 1}, 'alpha',sig_alpha, 'tail','left'); 
        if dfof_tg_avg(icell, iori, 1) < resp_threshold; sig_ori_driven(icell, iori) = 0; end
        
        ori_resp = [ori_resp; resp_cond{icell, iori, 1}];
        ori_list = [ori_list; iori * ones(length(resp_cond{icell, iori, 1}),1)];
        baseline = [baseline; base_cond{icell, iori, 1}];
    end
    
    ori_resp_list = [ori_resp, ori_list];
    baseline(:,2) = zeros(size(baseline,1),1);
    anova_compare = [ori_resp_list; baseline];
    data = anova_compare(:,1);
    label = cellstr(num2str(anova_compare(:,2)));
    
    [p, ~] = anova1(data,label,'off');
    if p < sig_alpha; sig_vis_driven(icell) = 1; end
end

end