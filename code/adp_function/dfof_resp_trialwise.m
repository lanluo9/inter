function [dfof_ad_trial, dfof_base_trial] = dfof_resp_trialwise(dfof_align, save_flag)

% input: dfof_align by ad. save_flag to save mat
% output: ncell x ntrial x (nstim + 1) response matrix for anova

global ncell nori nisi id_isi3 id_ori range_base range_resp

base_cond = cell(ncell, nori, nisi); resp_cond = cell(ncell, nori, nisi);
init = pi * ones(ncell, nori, nisi); dfof_avg = init; dfof_sem = init; dfof_std = init;

for icell = 1 : ncell
for iori = 1 : nori 
for iisi =  1 : length(id_isi3) 
    idx = intersect(id_ori{iori}, id_isi3{iisi}); 
    base_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_base)),2); 
    resp_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_resp)),2);

    ntrial_cond = length(base_cond{icell, iori, iisi});
    dfof_avg(icell, iori, iisi) = mean( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
    dfof_sem(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} ) / sqrt(ntrial_cond);
    dfof_std(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
end
end
end

dfof_tg_avg = dfof_avg;
dfof_tg_sem = dfof_sem;
dfof_tg_std = dfof_std;

if save_flag; save dfof_tg.mat dfof_tg_avg dfof_tg_sem dfof_tg_std; end

end