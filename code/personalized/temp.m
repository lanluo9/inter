t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_750,:)), 1));
t_ad_750 = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_250,:)), 1));
t_ad_250 = squeeze(nanmean(t(:,:), 1)); 

t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_750,:)), 1)); 
t_tg_750 = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_250,:)), 1)); 
t_tg_250 = squeeze(nanmean(t(:,:), 1)); 

plot(t_ad_750(1:range), 'r'); hold on; plot(t_ad_250(1:range), 'r--'); 
plot(t_tg_750(1:range), 'b'); plot(t_tg_250(1:range), 'g'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
legend('ad align 750', 'ad align 250', 'targ align 750', 'targ align 250')