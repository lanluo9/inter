by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);

clear adp_sig_area
for iarea = 1 : narea
    area_set_seq = by_area_id{iarea};
    adp_sig_area(iarea).dfof_equiv_ad_ns = [];
    adp_sig_area(iarea).dfof_equiv_ad = [];
    adp_sig_area(iarea).adp_vis_driven = [];
    adp_sig_area(iarea).adp_ns = [];
    
    for iset = 1 : length(area_set_seq)
        adp_sig_area(iarea).dfof_equiv_ad_ns = [adp_sig_area(iarea).dfof_equiv_ad_ns; ...
            adp_sig(area_set_seq(iset)).dfof_equiv_ad_ns];
        adp_sig_area(iarea).dfof_equiv_ad = [adp_sig_area(iarea).dfof_equiv_ad; ...
            adp_sig(area_set_seq(iset)).dfof_equiv_ad];
        
        adp_sig_area(iarea).adp_vis_driven = [adp_sig_area(iarea).adp_vis_driven; ...
            adp_sig(area_set_seq(iset)).adp_vis_driven(:,2)];
        adp_sig_area(iarea).adp_ns = [adp_sig_area(iarea).adp_ns; ...
            adp_sig(area_set_seq(iset)).adp_ns(:,2)];
    end
end

%%
for iarea = 1 : narea
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    histogram(adp_sig_area(iarea).dfof_equiv_ad, 'BinWidth',0.02); hold on; % nbin 25
    histogram(adp_sig_area(iarea).dfof_equiv_ad_ns, 'BinWidth',0.02)
    legend('vis driven resp ad', 'ns resp ad')
    
    subplot(2,1,2)
    histogram(adp_sig_area(iarea).adp_vis_driven, 'BinWidth',1); hold on; 
    histogram(adp_sig_area(iarea).adp_ns, 'BinWidth',1)  % nbin 50
    legend('vis driven adp', 'ns adp')
    
    saveas(gcf, ['ttest adp for area ', num2str(iarea)], 'jpg')
end
