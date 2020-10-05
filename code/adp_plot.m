%% set up 

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\mat

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806, ...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
for iset = 1 : length(dataset_list.date)
    areamouse_seq{iset} = [dataset_list.area{1,iset} '_' num2str(dataset_list.mouse(iset))];
end

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
set = cell(nset, 5);
for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset))
    imouse = ['i', mouse];
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];

    result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
    cd(result_folder{iset});

    mat_files = dir('*.mat'); 
    for imat = 1:length(mat_files) 
        set{iset, imat} = load(mat_files(imat).name); 
    end
end

%% adaptation index
% adp plots & trace: take only vis_driven_ad cells & only resp_targ0
% refer to Jin2020 Fig 3D & 4B

ndelta = 8;
ngap = 2;

adp = struct;
low_ad_resp = {};
for iset = 1 : nset
    cell_list_now = find(set{iset, 1}.vis_driven_ad); 
    
% with-adapter / no-adapter resp to same targ with same isi - not good bc adp is stim-specific
adp_ratio = zeros(length(cell_list_now), ndelta, ngap);
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    for idelta = 1:ndelta
    for igap = 1:ngap
        dfof_equiv_ad = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_noad = set{iset, 2}.dfof_avg_merge(icell, idelta, 1);
        adp_ratio(ii, idelta, igap) = dfof_equiv_ad / dfof_equiv_noad - 1;
    end
    end
end

% with-ad targ0 vs no-ad targ0
adp_ratio_targ0 = zeros(length(cell_list_now), 1, ngap);
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8; % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ(ii, igap) = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_noad_targ(ii, igap) = set{iset, 2}.dfof_avg_merge(icell, idelta, 1);
        adp_ratio_targ0(ii, 1, igap) = dfof_equiv_ad_targ(ii, igap) / dfof_equiv_noad_targ(ii, igap) - 1;
    end
end

% with-ad targ0 vs its own ad0
adp_ratio_a0t0 = zeros(length(cell_list_now), 1, ngap);
load(fullfile(result_folder{iset}, 'pre-processing', 'resp_ad.mat'))
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8; % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ(ii, igap) = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_ad(ii, igap) = dfof_avg_ad(icell, idelta);
        adp_ratio_a0t0(ii, 1, igap) = dfof_equiv_ad_targ(ii, igap) / dfof_equiv_ad(ii, igap) - 1;
    end
end

% histogram(dfof_equiv_ad(:,1),100); hold on; grid minor
% line([0.01, 0.01], [0, 4], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
% line([-0.01, -0.01], [0, 4], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
% xlabel('response to adapter')
% ylabel('number of cells')

low_ad_resp{iset} = logical(abs(dfof_equiv_ad(:,1))<0.01); % two col of dfof_equiv_ad are identical
% low_ad_resp{iset} = logical(abs(dfof_equiv_noad_targ(:,1))<0.01); 
nlow(iset) = sum(low_ad_resp{iset});

% adp(iset).adp_ratio = adp_ratio;
% adp(iset).adp_ratio_targ0 = adp_ratio_targ0;
% adp(iset).adp_ratio_a0t0 = adp_ratio_targ0;
adp(iset).adp_ratio_a0t0 = adp_ratio_a0t0;

end

for iset = 1 : nset
    adp_ratio_a0t0_high_ad = zeros(sum(~low_ad_resp{iset}),1,ngap);
    for igap = 1:ngap
        high_ad_resp = ~low_ad_resp{iset};
        adp_ratio_a0t0_high_ad(:,1,igap) = adp(iset).adp_ratio_a0t0(high_ad_resp, 1, igap);
    end
    adp(iset).adp_ratio_a0t0_high_ad = adp_ratio_a0t0_high_ad;
    
    for igap = 1:ngap
        adp(iset).adp_pos{igap} = logical(adp_ratio_a0t0_high_ad(:,igap)>0);
        adp(iset).adp_neg{igap} = logical(adp_ratio_a0t0_high_ad(:,igap)<0); % when adp inhibits response, as expected
    end
end

%% san check for a0t0

% t = adp(5).adp_ratio_a0t0(:,:,1);
% t2 = adp(5).adp_ratio_a0t0(:,:,2);
% histogram(t(abs(t)<40),100); hold on; histogram(-1*t2(abs(t2)<40),100)

% histogram(dfof_equiv_ad,100); grid on; grid minor
% histogram(dfof_equiv_ad_targ,100)
% 
% for iset = 1:nset
%     low_ad_resp(iset) = logical(sum(abs(dfof_equiv_ad(iset, :, :))<0.01,2)); 
%     nlow(iset) = sum(low_ad_resp)
%     dfof_equiv_ad(low_ad_resp)
%     dfof_equiv_ad_targ(low_ad_resp)
% end

%% adp by area
% violin plot and avg-sem
% use adp_ratio_a0t0 or adp_ratio_targ0? -> use a0t0

by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);

for iarea = 1 : narea
    area_set_seq = by_area_id{iarea};
    adp_ratio_now = cell(1,2);
    
    for igap = 1 : ngap
    for iset = 1 : length(area_set_seq)
        tt = adp(area_set_seq(iset)).adp_ratio_a0t0_high_ad(:,:,igap); 
        tt = mean(tt,2);
        tt = tt(:);
        adp_ratio_now{igap} = [adp_ratio_now{igap}; tt];
    end
    end
    adp_area(iarea).adp_ratio = adp_ratio_now;
    
    area_name = convertCharsToStrings(dataset_list.area{1,iarea});
    adp_area(iarea).name = repmat(area_name, length(adp_ratio_now{igap}), 1);
end

cd(result_prefix)
isi_str = {'isi 750', 'isi 250'};

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(1,2,igap);
    values = [adp_area(1).adp_ratio{1, igap}; adp_area(2).adp_ratio{1, igap};...
        adp_area(3).adp_ratio{1, igap}];
%     histogram(values,100); close
%     outlier = find( abs(values) - mean(values(:)) > 3*std(values(:)) );
%     values(outlier) = NaN;
    
    ncell_area = [length(adp_area(1).name); length(adp_area(2).name); length(adp_area(3).name)];

    names = [adp_area(1).name; adp_area(2).name; adp_area(3).name];
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
        tt = adp_area(iarea).adp_ratio{1, igap}; 
%         outlier = find( abs(tt) - mean(tt(:)) > 3*std(tt(:)) );
%         tt(outlier) = NaN;
        adp_area_avg(iarea, igap) = nanmean(tt);
        adp_area_ste(iarea, igap) = nanstd(tt) ./ sqrt(length(tt(~isnan(tt))));
    end
end

figure
color_list = {[0,0,1], [1,0,0]};
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

%% adp by mouse

by_mouse_id = {[1,2,3], [4,5,6], [7,8]}; nmouse = length(by_mouse_id);
adp_mouse = struct;

for imouse = 1 : nmouse
    mouse_set_seq = by_mouse_id{imouse};
    adp_ratio_now = cell(1,2);
    
    for igap = 1 : ngap
    for iset = 1 : length(mouse_set_seq)
        tt = adp(mouse_set_seq(iset)).adp_ratio_a0t0_high_ad(:,:,igap); 
        tt = mean(tt,2);
        tt = tt(:);
        adp_ratio_now{igap} = [adp_ratio_now{igap}; tt];
    end
    end
    adp_mouse(imouse).adp_ratio = adp_ratio_now;
    
    mouse_name = num2str(dataset_list.mouse(mouse_set_seq(iset))); mouse_name = convertCharsToStrings(mouse_name);
    adp_mouse(imouse).name = repmat(mouse_name, length(adp_ratio_now{igap}), 1);
end

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(1,2,igap);
    values = [adp_mouse(1).adp_ratio{1, igap}; adp_mouse(2).adp_ratio{1, igap};...
        adp_mouse(3).adp_ratio{1, igap}];
%     outlier = find( abs(values) - mean(values(:)) > 3*std(values(:)) );
%     values(outlier) = NaN;

    ncell_mouse = [length(adp_mouse(1).name); length(adp_mouse(2).name); length(adp_mouse(3).name)];
    names = [adp_mouse(1).name; adp_mouse(2).name; adp_mouse(3).name];
    
    vs = violinplot(values, names, 'GroupOrder', {'1322','1323','1324'}); hold on
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 3.5])
    line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
    for itext = 1 : length(ncell_mouse)
        text(itext, -3 + 0.5, ...
            ['n=', num2str(ncell_mouse(itext))], 'HorizontalAlignment', 'center')
    end
end
yl = ylim;
ylim(hAx,[-3, yl(2)])
% saveas(gcf, ['adp ratio a0t0 across mouse'], 'jpg'); close 

% ylim(hAx,[-10, 10])
% % saveas(gcf, ['adp ratio targ0 across mouse'], 'jpg'); close 

%%
for igap = 1:ngap
    for imouse = 1:nmouse
        tt = adp_mouse(imouse).adp_ratio{1, igap}; 
%         outlier = find( abs(tt) - mean(tt(:)) > 3*std(tt(:)) );
%         tt(outlier) = NaN;
        adp_mouse_avg(imouse, igap) = nanmean(tt);
        adp_mouse_ste(imouse, igap) = nanstd(tt) ./ sqrt(length(tt(~isnan(tt))));
    end
end

color_list = {[0,0,1], [1,0,0]};
for igap = 1:ngap
    h{igap,1} = scatter(1:3, adp_mouse_avg(:, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:3, adp_mouse_avg(:, igap), adp_mouse_ste(:, igap), 'color', color_list{igap}, 'LineStyle','none')
end
ylim([-1, 1])
yl = ylim;
for itext = 1 : length(ncell_mouse)
    text(itext, yl(1) + 0.05, ...
        ['n=', num2str(ncell_mouse(itext))], 'HorizontalAlignment', 'center')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:3]); xticklabels({'1322','1323','1324'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250', 'Location','northeast'); legend boxoff
xlim([0.5, 3.5])
% saveas(gcf, ['adp ratio a0t0 across mouse avg sem'], 'jpg'); close 
% % saveas(gcf, ['adp ratio targ0 across mouse avg sem'], 'jpg'); close 

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
        tt = adp(iset).adp_ratio_a0t0_high_ad(:,:,igap); 
        tt = mean(tt,2);
        adp_ratio_now{igap} = tt(:);
    end
    adp_set(iset).adp_ratio = adp_ratio_now;
    
    area_name = convertCharsToStrings(set_name{iset});
    adp_set(iset).name = repmat(area_name, length(adp_ratio_now{igap}), 1);
end

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(2,1,igap);
    values = []; ncell_set = []; names = [];
    
    for iset = 1:nset
        values = [values; adp_set(iset).adp_ratio{1, igap}];
%         outlier = find( abs(values) - mean(values(:)) > 3*std(values(:)) );
%         values(outlier) = NaN;
        
        ncell_set = [ncell_set; length(adp_set(iset).name)];
        names = [names; strrep(adp_set(iset).name, '_13', '.')];
    end
    
    vs = violinplot(values, names, 'GroupOrder', {'V1.22', 'V1.23', 'V1.24', ...
        'LM.22', 'LM.23', 'LM.24',  'LI.22', 'LI.23'}); hold on
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 8.5])
    line([0.5, 8.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
    
%     yl = ylim;
    for itext = 1 : length(ncell_set)
        text(itext, 6, ...
            ['n=', num2str(ncell_set(area_mouse_id(itext)))], 'HorizontalAlignment', 'center')
    end
end
ylim(hAx,[-3, 6.5])
% ylim(hAx,[-10, 10])
% saveas(gcf, ['adp ratio a0t0 across area-mouse zoom in'], 'jpg'); close 
% % saveas(gcf, ['adp ratio targ0 across area-mouse zoom in'], 'jpg'); close 

%%
for igap = 1:ngap
    for iset = 1:nset
        tt = adp_set(iset).adp_ratio{1, igap}; 
%         outlier = find( abs(tt) - mean(tt(:)) > 3*std(tt(:)) );
%         tt(outlier) = NaN;
        adp_set_avg(iset, igap) = nanmean(tt);
        adp_set_ste(iset, igap) = nanstd(tt) ./ sqrt(length(tt(~isnan(tt))));
    end
end

% color_list = {[24,95,173], [17,174,207], [176,111,175]}
color_list = {[0,0,1], [1,0,0]};
for igap = 1:ngap
    h{igap,1} = scatter(1:8, adp_set_avg(area_mouse_id, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:8, adp_set_avg(area_mouse_id, igap), adp_set_ste(area_mouse_id, igap), 'color', color_list{igap}, 'LineStyle','none')
end
ylim([-1.5, 1.5])
yl = ylim;
for itext = 1 : length(ncell_set)
    text(itext, yl(1) + 0.1, ...
        ['n=', num2str(ncell_set(area_mouse_id(itext)))], 'HorizontalAlignment', 'center')
end
line([0.5, 8.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:8]); xticklabels({'V1.1322', '23', '24', ...
        'LM.1322', '23', '24',  'LI.1322', '23'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250', 'Location','northwest'); legend boxoff
xlim([0.5, 8.5])
% saveas(gcf, ['adp ratio a0t0 across area-mouse avg sem'], 'jpg'); close
% % saveas(gcf, ['adp ratio targ0 across area-mouse avg sem'], 'jpg'); close

%% trace grand average for area or for mouse

trace = struct;
for iset = 1:nset
    trace(iset).ncell = size(set{iset, 5}.trace_targ0_750, 1) - nlow(iset);
    trace(iset).trace_avg{1,1} = nanmean(set{iset, 5}.trace_targ0_750(~low_ad_resp{iset},:),1);
    trace(iset).trace_avg{2,1} = nanmean(set{iset, 5}.trace_targ0_250(~low_ad_resp{iset},:),1);
%     trace(iset).trace_std{1,1} = nanstd(set{iset, 5}.trace_targ0_750(~low_ad_resp{iset},:),1);
%     trace(iset).trace_std{2,1} = nanstd(set{iset, 5}.trace_targ0_250(~low_ad_resp{iset},:),1);

%     trace(iset).ncell_neg = [sum(adp(iset).adp_neg{1}), sum(adp(iset).adp_neg{2})]; % why complimentary mirror?
%     trace(iset).ncell_pos = [sum(adp(iset).adp_pos{1}), sum(adp(iset).adp_pos{2})];
%     trace(iset).trace_avg{1} = nanmean(set{iset, 5}.trace_targ0_750(adp(iset).adp_pos{1},:),1);
%     trace(iset).trace_avg{2} = nanmean(set{iset, 5}.trace_targ0_250(adp(iset).adp_pos{2},:),1);
end

%% trace for area

by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);
figure('units','normalized','outerposition',[0 0 1 1]);

for igap = 1 : ngap
    subplot(1,2,igap)
    % igap = 1; % plot isi 750 -> 250
    trace_area_avg = []; trace_area_sets_std = [];

    for iarea = 1 : narea
        area_set_seq = by_area_id{iarea};
        trace_area_sets = []; 
        ncell_sum = 0;

        for iset = 1 : length(area_set_seq)
            trace_area_sets = [trace_area_sets; ...
                trace(area_set_seq(iset)).trace_avg{igap} .* trace(area_set_seq(iset)).ncell];
            ncell_sum = ncell_sum + trace(area_set_seq(iset)).ncell;
        end
        trace_area_avg(iarea,:) = sum(trace_area_sets,1) ./ ncell_sum;

    %     % possible bug: how to merge std?
    %     for iset = 1 : length(area_set_seq)
    %         trace_area_sets_std = [trace_area_sets_std; ...
    %             trace_neg(area_set_seq(iset)).trace_std{igap} ./ trace_neg(area_set_seq(iset)).ncell];
    %         ncell_sum = ncell_sum + trace_neg(area_set_seq(iset)).ncell;
    %     end
    %     trace_area_std(iarea,:) = sum(trace_area_sets_std,1); % .* ncell_sum;
    end

    trace_len = 3.5 * 30; % trace len 3.5 s according to Jin2019 Fig 2B
    color_list = {[0,0,1], [0,1,0], [1,0,0]};
    for iarea = 1 : narea
        plot(1:trace_len, trace_area_avg(iarea, 1:trace_len), 'color',color_list{iarea}); hold on
        xlim([0, 105])
    end
    
    % for iarea = 1 : narea
    %     x = 1:trace_len; x2 = [x, fliplr(x)];
    %     curve1 = trace_area_avg(iarea, 1:trace_len) + trace_area_std(iarea, 1:trace_len); curve2 = trace_area_avg(iarea, 1:trace_len) - trace_area_std(iarea, 1:trace_len);
    %     inBetween = [curve1, fliplr(curve2)];
    %     h = fill(x2, inBetween, color_list{iarea}, 'edgecolor','none'); 
    %     h.FaceAlpha = 0.3;
    % end
    
    legend('V1', 'LM', 'LI', 'Location','northeast'); legend boxoff
    title([isi_str{igap}])

end
% saveas(gcf, ['trace across area'], 'jpg'); close

%{
%% trace for mouse

by_mouse_id = {[1,2,3], [4,5,6], [7,8]}; nmouse = length(by_mouse_id);
igap = 2; 
trace_mouse_avg = []; trace_mouse_sets_std = [];

for imouse = 1 : nmouse
    mouse_set_seq = by_mouse_id{imouse};
    trace_mouse_sets = []; 
    ncell_sum = 0;
    
    for iset = 1 : length(mouse_set_seq)
        trace_mouse_sets = [trace_mouse_sets; ...
            trace(mouse_set_seq(iset)).trace_avg{igap} .* trace(mouse_set_seq(iset)).ncell];
        ncell_sum = ncell_sum + trace(mouse_set_seq(iset)).ncell;
    end
    trace_mouse_avg(imouse,:) = sum(trace_mouse_sets,1) ./ ncell_sum;
    
    % possible bug: how to merge std?
    for iset = 1 : length(mouse_set_seq)
        trace_mouse_sets_std = [trace_mouse_sets_std; ...
            trace(mouse_set_seq(iset)).trace_std{igap} ./ trace(mouse_set_seq(iset)).ncell];
        ncell_sum = ncell_sum + trace(mouse_set_seq(iset)).ncell;
    end
    trace_mouse_std(imouse,:) = sum(trace_mouse_sets_std,1); % .* ncell_sum;
    
end

trace_len = 3.5 * 30; % trace len 3.5 s according to Jin2019 Fig 2B
color_list = {[0,0,1], [0,1,0], [1,0,0]};
for imouse = 1 : nmouse
    plot(1:trace_len, trace_mouse_avg(imouse, 1:trace_len), 'color',color_list{imouse}); hold on
    xlim([0, 105])
end
legend('1322', '1323', '1324', 'Location','northeast'); legend boxoff
for imouse = 1 : nmouse
    x = 1:trace_len; x2 = [x, fliplr(x)];
    curve1 = trace_mouse_avg(imouse, 1:trace_len) + trace_mouse_std(imouse, 1:trace_len); curve2 = trace_mouse_avg(imouse, 1:trace_len) - trace_mouse_std(imouse, 1:trace_len);
    inBetween = [curve1, fliplr(curve2)];
    h = fill(x2, inBetween, color_list{imouse}, 'edgecolor','none'); 
    h.FaceAlpha = 0.3;
end
% % saveas(gcf, ['trace across mouse isi ', num2str(igap)], 'jpg'); close

for imouse = 1:3
    resp_ad(imouse) = max(trace_mouse_avg(imouse, 1:15)) - min(trace_mouse_avg(imouse, 1:15));
    resp_ad_targ(imouse) = max(trace_mouse_avg(imouse, 16:60)) - min(trace_mouse_avg(imouse, 16:40));
    adp_mouse_trace(imouse) = resp_ad_targ(imouse)/resp_ad(imouse) - 1;
end
adp_mouse_trace

%% trace for area but list mouse

by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);
igap = 1;

for iarea = 1 : narea
    area_set_seq = by_area_id{iarea};
    trace_area_sets = []; 
    ncell_sum = 0;
    
    for iset = 1 : length(area_set_seq)
        trace_area_sets = [trace_area_sets; ...
            trace(area_set_seq(iset)).trace_avg{igap,1}];
    end
    trace_area_mouse{iarea} = trace_area_sets;
    
end

area_str = {'V1', 'LM', 'LI'};
figure('units','normalized','outerposition',[0 0 1 1/2]);
for iarea = 1 : narea
    subplot(1,3,iarea)
    for imouse = 1:size(trace_area_mouse{1,iarea},1)
        plot(trace_area_mouse{1,iarea}(imouse, 1:trace_len)); hold on
    end
    xlim([0, 105])
    ylim([-0.05, 0.25])
    xticks(50); 
    xticklabels(area_str{iarea})
    if iarea == 3
        legend('1322', '1323', 'Location','northeast'); 
    else
        legend('1322', '1323', '1324', 'Location','northeast'); 
    end
    legend boxoff
end
% % % saveas(gcf, ['trace across area-mouse isi ', num2str(igap)], 'jpg'); close
%}

%% 3-way anova for area mouse isi & multiple comparison

cd C:\Users\lan\Documents\repos\inter\mat
adp_ratio_seq = []; ncell_set = [];
for igap = 1:ngap
    for iset = 1:nset
        append_seq = squeeze(adp(iset).adp_ratio_a0t0_high_ad(:,1,igap));
        ncell_set(iset) = length(append_seq);
        adp_ratio_seq = [adp_ratio_seq; append_seq];
    end
end

ncell_set2 = [ncell_set, ncell_set];
list_area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM',...
             'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
list_mouse = {'1322','1322','1322', '1323','1323','1323', '1324','1324',...
              '1322','1322','1322', '1323','1323','1323', '1324','1324'};
list_isi = {'750','750','750','750', '750','750','750','750',... 
            '250','250','250','250', '250','250','250','250'};
area_seq = char; mouse_seq = char; isi_seq = char; 
for iset = 1 : nset*2
    for icell = 1 : ncell_set2(iset)
        area_seq = [area_seq; list_area{iset}];
        mouse_seq = [mouse_seq; list_mouse{iset}];
        isi_seq = [isi_seq; list_isi{iset}];
    end
end

% p_anova_n = anovan(adp_ratio_seq,{area_seq mouse_seq isi_seq}, 'model',3, 'varnames',{'area', 'mouse', 'isi'})
% p_anova_n = anovan(adp_ratio_seq,{area_seq mouse_seq}, 'model',2, 'varnames',{'area', 'mouse'})
% p_anova_n = anovan(adp_ratio_seq,{area_seq isi_seq}, 'model',2, 'varnames',{'area', 'isi'})

[p,tbl,stats] = anovan(adp_ratio_seq,{area_seq mouse_seq isi_seq},'model','interaction',...
    'varnames',{'area', 'mouse', 'isi'});
% saveas(gcf, ['anovan area-mouse-isi'], 'jpg'); close

% results = multcompare(stats,'Dimension',[1 2])
% results = multcompare(stats,'Dimension',[1 3])
results = multcompare(stats,'Dimension',[1 2 3])
set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['multcompare area-mouse-isi'], 'png'); close

%% 
% take only vis_driven & good_fit cells
% refer to Jin2019 Fig 2D-G

