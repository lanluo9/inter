for icell = 61:83

t = squeeze(nanmean(squeeze(tc_trial_align_ad(vis_driven_cell_list(icell),:,:)), 1));
t_ad = squeeze(nanmean(t(:,:), 1)); 

plot(t_ad(1:50)); hold on; 

end

grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
% legend('ad align', 'targ align 750', 'targ align 250')