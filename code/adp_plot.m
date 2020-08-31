%% set up 

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\code

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806,...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
for iset = 1 : length(dataset_list.date)
    areamouse_seq{iset} = [dataset_list.area{1,iset} '_' num2str(dataset_list.mouse(iset))];
end

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
set = cell(nset, 4);
for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset))
    imouse = ['i', mouse];
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];

    result_prefix = 'C:\Users\lan\Documents\repos\inter\code\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
    cd(result_folder{iset});

    mat_files = dir('*.mat'); 
    for imat = 1:length(mat_files) 
        set{iset, imat} = load(mat_files(imat).name); 
    end
end

%% adaptation index

ndelta = 8;
ngap = 2;

adp = struct;
for iset = 1 : nset
    
% with-adapter / no-adapter resp to same targ with same isi
cell_list_now = find(set{iset, 1}.vis_driven); 
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
% cell_list_now = set{iset, 1}.pref_0_cell; % only a subpop. perhaps doesn't matter bc Jin2019 Fig1D
adp_ratio_targ0 = zeros(length(cell_list_now), ngap);

for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8; % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_noad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, 1);
%         dfof_equiv_noad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 9:11)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 1:3)),2);
        adp_ratio_targ0(ii, igap) = dfof_equiv_ad_targ / dfof_equiv_noad_targ - 1;
    end
end

% with-ad targ0 vs its own ad0
adp_ratio_a0t0 = zeros(length(cell_list_now), ngap);
load(fullfile(result_folder{iset}, 'pre-processing', 'resp_ad.mat'))
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    idelta = 8 % targ0 only! adp is ori-specific
    for igap = 1:ngap
        dfof_equiv_ad_targ = set{iset, 2}.dfof_avg_merge(icell, idelta, igap+1); %750-250
        dfof_equiv_ad = dfof_avg_ad(icell, idelta);
        adp_ratio_a0t0(ii, igap) = dfof_equiv_ad_targ / dfof_equiv_ad - 1;
    end
end

adp(iset).adp_ratio = adp_ratio;
adp(iset).adp_ratio_targ0 = adp_ratio_targ0;
adp(iset).adp_ratio_a0t0 = adp_ratio_a0t0;

end

%% adp by area
% violin plot and avg-sem

by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);

for iarea = 1 : narea
    area_set_seq = by_area_id{iarea};
    adp_ratio_now = cell(1,2);
    
    for igap = 1 : ngap
    for iset = 1 : length(area_set_seq)
        tt = adp(area_set_seq(iset)).adp_ratio(:,:,igap); 
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
isi_str = {'isi 750', 'isi 250'}

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(1,2,igap)
    values = [adp_area(1).adp_ratio{1, igap}; adp_area(2).adp_ratio{1, igap};...
        adp_area(3).adp_ratio{1, igap}];
    outlier = find(values > mean(values(:) + std(values(:))) ...
        | values < mean(values(:) - std(values(:))));
    values(outlier) = NaN;

    names = [adp_area(1).name; adp_area(2).name; adp_area(3).name];
    vs = violinplot(values, names, 'GroupOrder', {'V1','LM','LI'});
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 3.5])
    line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
end
% saveas(gcf, ['adp idx across mouse w outlier'], 'jpg'); 
% ylim(hAx,[-4, 4])
% saveas(gcf, ['adp idx across mouse no outlier'], 'jpg'); 
% ylim(hAx,[-1, 1])
% saveas(gcf, ['adp idx across mouse zoom in'], 'jpg'); close 

for igap = 1:ngap
    for iarea = 1:narea
        tt = adp_area(iarea).adp_ratio{1, igap}; 
        outlier = find(tt > mean(tt(:) + std(tt(:))) | tt < mean(tt(:) - std(tt(:))));
        tt(outlier) = NaN;
        adp_area_avg(iarea, igap) = nanmean(tt);
        adp_area_ste(iarea, igap) = nanstd(tt) ./ length(tt(~isnan(tt)));
    end
end

% color_list = {[24,95,173], [17,174,207], [176,111,175]}
color_list = {[0,0,1], [1,0,0]};
for igap = 1:ngap
    h{igap,1} = scatter(1:3, adp_area_avg(:, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:3, adp_area_avg(:, igap), adp_area_ste(:, igap), 'color', color_list{igap}, 'LineStyle','none')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:3]); xticklabels({'V1', 'LM', 'LI'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250'); legend boxoff
xlim([0.5, 3.5])
ylim([-0.6, 0.1])
% saveas(gcf, ['adp idx across area avg sem'], 'jpg'); close 

%% adp by mouse

by_mouse_id = {[1,2,3], [4,5,6], [7,8]}; nmouse = length(by_mouse_id);
adp_mouse = struct;

for imouse = 1 : nmouse
    mouse_set_seq = by_mouse_id{imouse};
    adp_ratio_now = cell(1,2);
    
    for igap = 1 : ngap
    for iset = 1 : length(mouse_set_seq)
        tt = adp(mouse_set_seq(iset)).adp_ratio(:,:,igap); 
        tt = mean(tt,2);
        tt = tt(:);
        adp_ratio_now{igap} = [adp_ratio_now{igap}; tt];
    end
    end
    adp_mouse(imouse).adp_ratio = adp_ratio_now;
    
    mouse_name = num2str(dataset_list.mouse(mouse_set_seq(iset))); mouse_name = convertCharsToStrings(mouse_name);
    adp_mouse(imouse).name = repmat(mouse_name, length(adp_ratio_now{igap}), 1);
end

cd(result_prefix)
isi_str = {'isi 750', 'isi 250'}

figure('units','normalized','outerposition',[0 0 1 1]);
for igap = 1:ngap
    hAx(igap) = subplot(1,2,igap);
    values = [adp_mouse(1).adp_ratio{1, igap}; adp_mouse(2).adp_ratio{1, igap};...
        adp_mouse(3).adp_ratio{1, igap}];
    outlier = find(values > mean(values(:) + std(values(:))) ...
        | values < mean(values(:) - std(values(:))));
    values(outlier) = NaN;

    names = [adp_mouse(1).name; adp_mouse(2).name; adp_mouse(3).name];
    
    vs = violinplot(values, names, 'GroupOrder', {'1322','1323','1324'});
    ylabel(['adaptation index: ', isi_str{igap}]);
    xlim([0.5, 3.5])
    line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
end
% saveas(gcf, ['adp idx across mouse w outlier'], 'jpg'); 
% ylim(hAx,[-4, 4])
% saveas(gcf, ['adp idx across mouse no outlier'], 'jpg'); 
% ylim(hAx,[-1, 1])
% saveas(gcf, ['adp idx across mouse zoom in'], 'jpg'); close 

for igap = 1:ngap
    for imouse = 1:nmouse
        tt = adp_mouse(imouse).adp_ratio{1, igap}; 
        outlier = find(tt > mean(tt(:) + std(tt(:))) | tt < mean(tt(:) - std(tt(:))));
        tt(outlier) = NaN;
        adp_mouse_avg(imouse, igap) = nanmean(tt);
        adp_mouse_ste(imouse, igap) = nanstd(tt) ./ length(tt(~isnan(tt)));
    end
end

% color_list = {[24,95,173], [17,174,207], [176,111,175]}
color_list = {[0,0,1], [1,0,0]};
for igap = 1:ngap
    h{igap,1} = scatter(1:3, adp_mouse_avg(:, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:3, adp_mouse_avg(:, igap), adp_mouse_ste(:, igap), 'color', color_list{igap}, 'LineStyle','none')
end
line([0.5, 3.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:3]); xticklabels({'1322','1323','1324'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250'); legend boxoff
xlim([0.5, 3.5])
ylim([-0.7, 0.1])
% saveas(gcf, ['adp idx across mouse avg sem'], 'jpg'); close 

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
    h{igap,1} = scatter(1:8, adp_set_avg(area_mouse_id, igap), 5, color_list{igap}, 'filled'); hold on
    errorbar(1:8, adp_set_avg(area_mouse_id, igap), adp_set_ste(area_mouse_id, igap), 'color', color_list{igap}, 'LineStyle','none')
end
line([0.5, 8.5], [0, 0], 'Color', [0.7,0.7,0.7], 'LineWidth',1, 'LineStyle','--');
xticks([1:8]); xticklabels({'V1.1322', '23', '24', ...
        'LM.1322', '23', '24',  'LI.1322', '23'})
ylabel(['adaptation index']);
legend([h{1,1},h{2,1}], 'isi 750', 'isi 250'); legend boxoff
xlim([0.5, 8.5])
% saveas(gcf, ['adp idx across area-mouse avg sem'], 'jpg'); close

%% trace grand average for area or for mouse

