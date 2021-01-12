%%
area_merge = []; ori_pref = []; ori_orth = []; 
dfof_pref = []; dfof_orth = [];
ori_round_to = 0 : 22.5 : 157.5;

for iset = 1 : nset
    well_fit_cell = set(iset).cell_property.vis_cell_noad_tg & set(iset).cell_property.well_fit_cell;
    ncell_set(iset) = sum(well_fit_cell); 
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set(iset),1) * areacode;
    area_merge = [area_merge; temp];
    
    temp = set(iset).cell_property.ori_pref(well_fit_cell, 1);
    temp = interp1(ori_round_to,ori_round_to,temp, 'nearest','extrap');
    ori_pref = [ori_pref; temp];
    temp2 = temp + 90; temp2(temp2>=180) = temp2(temp2>=180) - 180;
    ori_orth = [ori_orth; temp2];
    
    for icell = 1 : ncell_set(iset)
        cell_idx = max(find(well_fit_cell, icell));
        ori_pref_idx = find(ori_round_to == temp(icell));
        ori_orth_idx = find(ori_round_to == temp2(icell));
        dfof_pref_cell = set(iset).dfof.dfof_tg(cell_idx, ori_pref_idx, 1);
        dfof_orth_cell = set(iset).dfof.dfof_tg(cell_idx, ori_orth_idx, 1);
        dfof_pref = [dfof_pref; dfof_pref_cell];
        dfof_orth = [dfof_orth; dfof_orth_cell];
    end
end
OSI = (dfof_pref - dfof_orth) ./ (dfof_pref + dfof_orth);
well_fit_flag = OSI; well_fit_flag(:) = 1;

%%

area_merge_not = []; ori_pref = []; ori_orth = []; 
dfof_pref = []; dfof_orth = [];

for iset = 1 : nset
    not_well_fit_cell = set(iset).cell_property.vis_cell_noad_tg & ~set(iset).cell_property.well_fit_cell;
    ncell_set(iset) = sum(not_well_fit_cell); 
    areacode = dataset_list.areacode{iset};
    temp = ones(ncell_set(iset),1) * areacode;
    area_merge = [area_merge; temp];
    
    temp = set(iset).cell_property.ori_pref(not_well_fit_cell, 1);
    temp = interp1(ori_round_to,ori_round_to,temp, 'nearest','extrap');
    ori_pref = [ori_pref; temp];
    temp2 = temp + 90; temp2(temp2>=180) = temp2(temp2>=180) - 180;
    ori_orth = [ori_orth; temp2];
    
    for icell = 1 : ncell_set(iset)
        cell_idx = max(find(not_well_fit_cell, icell));
        ori_pref_idx = find(ori_round_to == temp(icell));
        ori_orth_idx = find(ori_round_to == temp2(icell));
        dfof_pref_cell = set(iset).dfof.dfof_tg(cell_idx, ori_pref_idx, 1);
        dfof_orth_cell = set(iset).dfof.dfof_tg(cell_idx, ori_orth_idx, 1);
        dfof_pref = [dfof_pref; dfof_pref_cell];
        dfof_orth = [dfof_orth; dfof_orth_cell];
    end
end
OSI_not = (dfof_pref - dfof_orth) ./ (dfof_pref + dfof_orth);
well_fit_flag_not = OSI_not; well_fit_flag_not(:) = 0;

area = [area_merge; area_merge_not]; well_fit = [well_fit_flag; well_fit_flag_not];
OSI_all = [OSI; OSI_not];

% save OSI_area_vis.mat area well_fit OSI_all


% save set.mat set

