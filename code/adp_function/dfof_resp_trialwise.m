function [resp_ad, resp_base] = dfof_resp_trialwise(dfof_align, save_flag)

% input: dfof_align by ad or tg, linked w/ resp_win 'ad' or 'tg'. save_flag to save mat
% output: ncell x ntrial x (nstim + 1) response matrix for anova

global ncell nori nisi id_isi3 id_ori id_ad range_base range_resp

base = cell(ncell, 1); resp = cell(ncell, 1);
init = pi * ones(ncell, 1); dfof_avg = init; dfof_sem = init; dfof_std = init;

for icell = 1 : ncell
    base{icell} = mean(squeeze(dfof_align(icell, id_ad, range_base)),2); 
    resp{icell} = mean(squeeze(dfof_align(icell, id_ad, range_resp)),2);

    ntrial_ad = length(base{icell});
    dfof_avg(icell) = mean( resp{icell} - base{icell} );
    dfof_sem(icell) = std( resp{icell} - base{icell} ) / sqrt(ntrial_ad);
    dfof_std(icell) = std( resp{icell} - base{icell} );
end

dfof_ad_avg = dfof_avg;
dfof_ad_sem = dfof_sem;
dfof_ad_std = dfof_std;

if save_flag; save dfof_ad.mat dfof_ad_avg dfof_ad_sem dfof_ad_std; end

end