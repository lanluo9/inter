function dfof_resp(dfof_align, resp_win)

% 

global ncell nori nisi id_isi ori_seq ori_list id_ad id_noad_isi

switch resp_win
case 'ad'
    
case 'tg'
    
base_cond = cell(ncell, nori, nisi); resp_cond = cell(ncell, nori, nisi);
for icell = 1 : ncell
for iori = 1 : nori 
    id_ori = find(ori_seq == ori_list(iori));    
for iisi =  1 : length(id_noad_isi) 
    idx = intersect(intersect(id_noad_isi{iisi}, id_ori), id_ad); 
    base_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_base)),2); 
    resp_cond{icell, iori, iisi} = mean(squeeze(dfof_align(icell, idx, range_resp)),2);
end
end
end
base_avg_cond = pi * ones(ncell, nori, nisi); resp_avg_cond = pi * ones(ncell, nori, nisi); 
dfof_avg_cond = pi * ones(ncell, nori, nisi); dfof_ste_cond = pi * ones(ncell, nori, nisi); 
for icell = 1 : ncell
for iori = 1 : nori 
for iisi =  1 : nisi
    ntrial_cond = length(base_cond{icell, iori, iisi});
    base_avg_cond(icell, iori, iisi) = mean(base_cond{icell, iori, iisi}); 
    resp_avg_cond(icell, iori, iisi) = mean(resp_cond{icell, iori, iisi});
    dfof_avg_cond(icell, iori, iisi) = mean( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} );
    dfof_ste_cond(icell, iori, iisi) = std( resp_cond{icell, iori, iisi} - base_cond{icell, iori, iisi} ) / sqrt(ntrial_cond);
end
end
end
% save dfof_ad_tg.mat dfof_avg_cond dfof_ste_cond


dfof_avg_merge = cat(3, dfof_avg_noad, dfof_avg_750, dfof_avg_250);
dfof_ste_merge = cat(3, dfof_ste_noad, dfof_ste_750, dfof_ste_250);

cd(result_folder)
% save dfof_noad_750_250.mat dfof_avg_merge dfof_ste_merge

end