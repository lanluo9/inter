%% ori tuning & fit von Mises function for single cell (750)

dfof_avg_merge = cat(3, dfof_avg, dfof_avg_750, dfof_avg_250);
dfof_ste_merge = cat(3, dfof_ste, dfof_ste_750, dfof_ste_250);
fit_param_merge = cat(3, fit_param, fit_param_750, fit_param_250);

theta_finer = deg2rad(0:1:179);
subplot_title = {'control', 'isi 750 ms', 'isi 250 ms'};
vis_driven_ad_cell_list = find(vis_driven_ad); % using adapter-vis-driven cell bc more numerous

for iviscell = 1 : length(vis_driven_ad_cell_list)
    icell = vis_driven_ad_cell_list(iviscell);
    figure('units','normalized','outerposition',[0 0 1/2 1]);
for row = 1 : 3
    subplot(3,1,row)
    
    dfof_avg_now = dfof_avg_merge(:,:,row);
    dfof_ste_now = dfof_ste_merge(:,:,row);
    fit_param_now = fit_param_merge(:,:,row);
    
    errorbar([0,delta_list], [dfof_avg_now(icell,end), dfof_avg_now(icell,:)], ...
        [dfof_ste_now(icell,end), dfof_ste_now(icell,:)], 'LineStyle','none')
    hold on
    scatter([0,delta_list], [dfof_avg_now(icell,end), dfof_avg_now(icell,:)], 'b')

    t = num2cell(fit_param_now(icell, 2:end)); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});    
    y_fit(row,:) = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
    plot(rad2deg(theta_finer), y_fit(row,:), 'LineWidth', 1)
    
    ori_pref = rad2deg(u1_hat);
    ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref > 180) = ori_pref(ori_pref > 180) - 180;
    scatter(ori_pref, b_hat + R1_hat, 'r*') % mark pref ori of fit
    
    xlim([0-10, 180+10])
    tempmin = dfof_avg_merge - dfof_ste_merge; tempmin = tempmin(icell, :, :); ymin = min(tempmin(:));
    tempmax = dfof_avg_merge + dfof_ste_merge; tempmax = tempmax(icell, :, :); ymax = max(tempmax(:));
    padding = (ymax - ymin) ./ 50;
    ylim([ymin - padding, ymax + padding])
    
    xlabel('orientation (deg)')
    ylabel('dF/F')
    title(subplot_title{row})
end
    saveas(gcf, ['ori tuning across cond cell ', num2str(icell)], 'jpg')
    close
end


%% bootstrap -> goodness of fit 
% 
% % load ori_across_cells.mat
% nrun = 1000;
% theta = deg2rad(delta_list);
% 
% dfof_avg_runs = pi * ones(ncell, ndelta, nrun);
% dfof_ste_runs = pi * ones(ncell, ndelta, nrun);
% fit_param_runs = pi * ones(ncell, 7, nrun);
% ori_pref_runs = pi * ones(ncell, nrun);
% 
% for irun = 1 : nrun
%     irun
% 
%     for icell = 1 : ncell
%         
%         for idelta = 1 : ndelta
%             temp_win = cp_win_750{icell, idelta};
%             ntrial_cond = size(temp_win,1);
%             idx = 1 : ntrial_cond;
%             bootstrap_draw = round(ntrial_cond * 0.7);
%             idx_run = randsample(idx, bootstrap_draw, 1); % w replacement
% 
%             base_win = temp_win(idx_run,1); 
%             resp_win = temp_win(idx_run,2);
%             dfof_avg_runs(icell, idelta, irun) = mean( (resp_win - base_win) ./ mean(base_win) );
%             dfof_ste_runs(icell, idelta, irun) = std( (resp_win - base_win) ./ mean(base_win) ) ./ sqrt(ntrial_cond);
%         end
%         
%         data = dfof_avg_runs(icell, :, irun); 
%         [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
%         fit_param_runs(icell, :, irun) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%     %   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2
% 
%         ori_pref = rad2deg(u1_hat);
%         ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
%         ori_pref(ori_pref > 180) = ori_pref(ori_pref > 180) - 180;
%         ori_pref_runs(icell, irun) = ori_pref;
%     end
% end
% 
% % save ori_across_bootstrap_runs_with_replace.mat dfof_avg_runs dfof_ste_runs fit_param_runs ori_pref_runs
% 
% % sanity check
% tt = mean(dfof_avg_runs, 3);
% subplot(1,2,1)
% imagesc(dfof_avg_750); colorbar
% subplot(1,2,2)
% imagesc(mean(dfof_avg_runs, 3)); colorbar
% set(gcf, 'Position', get(0, 'Screensize'));
% 
