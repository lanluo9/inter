%% response amplitude by area
% 0. only vis_driven (by ad or tg) cells are included
% 1. resp_ad dist by area
% 2. resp_tg (averaged over oris) dist by area
% 3. radar plot of avg sem resp_tg w aligned max_resp_ori by area

area_ad = []; resp_ad = []; area_tg = []; resp_tg = [];

for iset = 1 : nset
    vis_ad_cell = set(iset).cell_property.vis_cell_ad; ncell_ad_set(iset) = sum(vis_ad_cell); 
    areacode = dataset_list.areacode{iset}; tt = ones(ncell_ad_set(iset),1) * areacode;
    area_ad = [area_ad; tt];
    tt = set(iset).dfof.dfof_ad(vis_ad_cell, 1);
    resp_ad = [resp_ad; tt];
    
    vis_tg_cell = set(iset).cell_property.vis_cell_noad_tg; ncell_tg_set(iset) = sum(vis_tg_cell); 
    areacode = dataset_list.areacode{iset}; tt = ones(ncell_tg_set(iset),1) * areacode;
    area_tg = [area_tg; tt];
    tt = set(iset).dfof.dfof_tg(vis_tg_cell, :, 1);
    resp_tg = [resp_tg; tt];
end

resp_tg_collapse_ori = mean(resp_tg, 2);
save amp_area.mat area_ad resp_ad area_tg resp_tg resp_tg_collapse_ori