by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);
by_dis_id = {[8], [1,7], [2,6], [3,5], [4]}; ndis = length(by_dis_id);

metrics_area = struct;
for iarea = 1 : narea
    area_set_seq = by_area_id{iarea};
    metrics_area_now = cell(1,2);
    
    for igap = 1 : ngap
    for iset = 1 : length(area_set_seq)
        tt = diff_dis(area_set_seq(iset)).resp_diff_dis(:,:,igap); 
        tt(abs(tt)>=1) = NaN;
        
        for idis = 1 : ndis
            dis_set_seq = by_dis_id{idis};
            
            for idelta = 1 : length(dis_set_seq)
                tt_seq = cat(1, [], tt(:,dis_set_seq(idelta)));
                tt_median(idis) = nanmedian(tt_seq, 1);
                tt_avg(idis) = nanmean(tt_seq, 1);
            end
        end
        metrics_area_now{igap} = tt_median;
    end
    end
    metrics_area(iarea).resp_diff_dis = metrics_area_now;
end

%%
dis_axis = sort(unique(dis_list));
figure('units','normalized','outerposition',[0 0 1 1]);
for iarea = 1:narea
    hAx(iarea) = subplot(1,narea,iarea);
    t_area = metrics_area(iarea).resp_diff_dis;
    for igap = 1:ngap
        plot(dis_axis, t_area{1,igap}); hold on
        yline(0)
    end
end
ylim(hAx,[-0.7, 0.21])
saveas(gcf, ['resp diff dis across areas median'], 'jpg'); close 

%%
cd(result_prefix)
isi_str = {'isi 750', 'isi 250'};

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(1,2,igap);
    values = [metrics_area(1).adp_ratio{1, igap}; metrics_area(2).adp_ratio{1, igap};...
        metrics_area(3).adp_ratio{1, igap}];
%     histogram(values,100); close
%     outlier = find( abs(values) - mean(values(:)) > 3*std(values(:)) );
%     values(outlier) = NaN;
    
    ncell_area = [length(metrics_area(1).name); length(metrics_area(2).name); length(metrics_area(3).name)];

    names = [metrics_area(1).name; metrics_area(2).name; metrics_area(3).name];
    vs = violinplot(values, names, 'GroupOrder', {'V1','LM','LI'}); hold on
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 3.5])
    line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
    
%     yl = ylim;
    for itext = 1 : length(ncell_area)
        text(itext, -3 + 0.5, ...
            ['n=', num2str(ncell_area(itext))], 'HorizontalAlignment', 'center')
    end
end
% ylim(hAx,[-2, 2])
% saveas(gcf, ['adp ratio a0t0 across areas violin'], 'jpg'); close 

%%
for igap = 1:ngap
    for iarea = 1:narea
        tt = metrics_area(iarea).adp_ratio{1, igap}; 
%         outlier = find( abs(tt) - mean(tt(:)) > 3*std(tt(:)) );
%         tt(outlier) = NaN;
        adp_area_avg(iarea, igap) = nanmean(tt);
        adp_area_ste(iarea, igap) = nanstd(tt) ./ sqrt(length(tt(~isnan(tt))));
    end
end

figure
color_list = {[0,0,1], [1,0,0]};
clear h
for igap = 1:ngap
    h{igap,1} = scatter(1:3, adp_area_avg(:, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:3, adp_area_avg(:, igap), adp_area_ste(:, igap), 'color', color_list{igap}, 'LineStyle','none')
end
ylim([-1.5,1.5])
yl = ylim;
for itext = 1 : length(ncell_area)
    text(itext, yl(1) + 0.1, ...
        ['n=', num2str(ncell_area(itext))], 'HorizontalAlignment', 'center')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:3]); xticklabels({'V1', 'LM', 'LI'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250', 'Location','northeast'); legend boxoff
xlim([0.5, 3.5])
% saveas(gcf, ['adp ratio a0t0 across areas avg sem'], 'jpg'); close 
% % saveas(gcf, ['adp ratio targ0 across areas avg sem'], 'jpg'); close 


%%

y = meshgrid(1:5);
rng default; % For reproducibility
y = y + normrnd(0,1,5,5);
p = anova1(y)

