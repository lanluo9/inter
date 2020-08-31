%% adp by area but list mouse

area_mouse_id = [1,4,7, 2,5,8, 3,6]; % set id sorted by area then by mouse
set_name = cell(nset, 1);
for iset = 1 : nset
    mouse = num2str(dataset_list.mouse(iset)); area = dataset_list.area{1,iset};
    areamouse = [area '_' mouse];
    set_name{iset} = areamouse;
end

for i = 1 : length(area_mouse_id)
    iset = area_mouse_id(i);
    adp_ratio_now = cell(1,2);
    
    for igap = 1 : ngap
        tt = adp(iset).adp_ratio(:,:,igap); 
        tt = mean(tt,2);
        adp_ratio_now{igap} = tt(:);
    end
    adp_set(iset).adp_ratio = adp_ratio_now;
    
    area_name = convertCharsToStrings(set_name{iset});
    adp_set(iset).name = repmat(area_name, length(adp_ratio_now{igap}), 1);
end


cd(result_prefix)
isi_str = {'isi 750', 'isi 250'};

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(2,1,igap);
    values = []; names = [];
    
    for iset = 1:nset
        values = [values; adp_set(iset).adp_ratio{1, igap}];
        outlier = find(values > 10 | values < -10);
%         outlier = find(values > mean(values(:) + std(values(:))) | values < mean(values(:) - std(values(:))));
        values(outlier) = NaN;
        names = [names; adp_set(iset).name];
    end
    
    vs = violinplot(values, names, 'GroupOrder', {'V1_1322', 'V1_1323', 'V1_1324', ...
        'LM_1322', 'LM_1323', 'LM_1324',  'LI_1322', 'LI_1323'});
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 8.5])
    line([0.5, 8.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
end
ylim(hAx,[-10, 10])
% saveas(gcf, ['adp idx across area-mouse w outlier zoom in'], 'jpg'); 
% ylim(hAx,[-5, 5])
% saveas(gcf, ['adp idx across area-mouse no outlier'], 'jpg'); 
% ylim(hAx,[-2, 2])
% saveas(gcf, ['adp idx across area-mouse no outlier zoom in'], 'jpg'); close 

%%
for igap = 1:ngap
    for iset = 1:nset
        tt = adp_set(iset).adp_ratio{1, igap}; 
%         outlier = find(tt > mean(tt(:) + std(tt(:))) | tt < mean(tt(:) - std(tt(:))));
        outlier = find(tt > 10 | tt < -10);
        tt(outlier) = NaN;
        adp_set_avg(iset, igap) = nanmean(tt);
        adp_set_ste(iset, igap) = nanstd(tt) ./ length(tt(~isnan(tt)));
    end
end

% color_list = {[24,95,173], [17,174,207], [176,111,175]}
color_list = {[0,0,1], [1,0,0]};
for igap = 1:ngap
    h{igap,1} = scatter(1:3, adp_set_avg(:, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:3, adp_set_avg(:, igap), adp_set_ste(:, igap), 'color', color_list{igap}, 'LineStyle','none')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:3]); xticklabels({'V1', 'LM', 'LI'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250'); legend boxoff
xlim([0.5, 3.5])
ylim([-0.6, 0.1])
% saveas(gcf, ['adp idx across area-mouse avg sem'], 'jpg'); close