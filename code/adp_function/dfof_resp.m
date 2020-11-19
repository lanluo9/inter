function [dfof_avg_cond, dfof_ste_cond] = dfof_resp(dfof_align, resp_win, save_flag)

% input: dfof_align by ad or tg, linked w/ resp_win 'ad' or 'tg'. save_flag to save mat
% output: dfof avg & sem 
% for ad = ncell x 1
% for tg = ncell x nori x nisi3 [noad 750 250]

global ncell nori nisi id_isi3 id_ori id_ad range_base range_resp

switch resp_win
case 'ad'
    base = cell(ncell, 1); resp = cell(ncell, 1);
    dfof_avg_cond = pi * ones(ncell, 1); dfof_ste_cond = pi * ones(ncell, 1); 

    for icell = 1 : ncell
        base{icell} = mean(squeeze(dfof_align(icell, id_ad, range_base)),2); 
        resp{icell} = mean(squeeze(dfof_align(icell, id_ad, range_resp)),2);

        ntrial_ad = length(base{icell});
        dfof_avg_cond(icell) = mean( resp{icell} - base{icell} );
        dfof_ste_cond(icell) = std( resp{icell} - base{icell} ) / sqrt(ntrial_ad);
    end

    if save_flag; save dfof_ad.mat dfof_avg_cond dfof_ste_cond; end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'tg'
    base_cond = cell(ncell, nori, nisi); resp_cond = cell(ncell, nori, nisi);
    dfof_avg_cond = pi * ones(ncell, nori, nisi); dfof_ste_cond = pi * ones(ncell, nori, nisi); 

    for icell = 1 : ncell
    for iori = 1 : nori 
    for iisi =  1 : length(id_isi3) 
        idx = intersect(id_ori{iori}, id_isi3{iisi}); 
        base_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_base)),2); 
        resp_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_resp)),2);

        ntrial_cond = length(base_cond{icell, iori, iisi});
        dfof_avg_cond(icell, iori, iisi) = mean( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
        dfof_ste_cond(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} ) / sqrt(ntrial_cond);
    end
    end
    end

    if save_flag; save dfof_tg.mat dfof_avg_cond dfof_ste_cond; end

end