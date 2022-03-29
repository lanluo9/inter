ccc;
area_names = {'LM','AL','PM','AM'};
summary_table = readtable('Z:\home\jen\Notebook files\HVA recording\Interneuron_properties_ratio.xlsx');

is_vClamp = cell2mat(cellfun(@(x) ~isempty(x),summary_table.HFS,'un',0));
filenames = cellfun(@(x,y) [x,y],summary_table.Date,summary_table.HFS,'un',0);

try 
    do_analysis = cell2mat(cellfun(@isempty,summary_table.Analyzed_ratio,'un',0));
catch
    do_analysis = isnan(summary_table.Analyzed_ratio);
end

keep_files = logical(do_analysis.*is_vClamp);
if any(~do_analysis)
    temp_disp = filenames(~do_analysis);
    disp(temp_disp(:));
    disp('found analysis files')
end


cellfun(@(x) int_pyr_ratio(x,false),filenames(keep_files),'un',0);

%%
clearvars -except PV_* SOM_*; clc;
area_names = {'LM','AL','PM','AM'};
setMin = true;
earlyLate_bound = 150;
summary_table = readtable('Z:\home\jen\Notebook files\HVA recording\Interneuron_properties_ratio.xlsx');

cell_name = 'Pyr';
if strcmp(cell_name,'PV')
    int_color = [0.3 0.5 0.8];
else
    int_color = [0.3 .8 0.4];
end
cmap = brewermap(4,'Set2')+0.01;
cmap_cell = mat2cell(cmap,[1 1 1 1],3);
keep_i = find(cell2mat(cellfun(@(x) contains(x,'X'),summary_table.KeepRatio,'un',0)));
type_idx = find(sum(cell2mat(cellfun(@(x) contains(x,cell_name),[summary_table.HS1 summary_table.HS2],'un',0)),2));
keep_i = intersect(keep_i,type_idx);

EPSC_grouped = cell(4,1);
param_grouped =  cell(4,1);
collapse_minCurrent = [];
for area_i = 1:4
    idx = find(cell2mat(cellfun(@(x) contains(x,area_names{area_i}),summary_table.Area,'un',0)));
    idx = intersect(idx,keep_i);
    
   
    for file_i = 1:size(idx,1)
        load([summary_table.Save_path_ratio{idx(file_i)},'.mat']);
        EPSC_grouped{area_i} = [EPSC_grouped{area_i} EPSC];
        param_grouped{area_i} = [param_grouped{area_i} param];
    end
    time_toPeak{area_i} = vertcat(EPSC_grouped{area_i}.timeToPeak)/20;
    int_pyr_ratio_grouped{area_i} = [EPSC_grouped{area_i}.Int_Pyr_ratio];
    filenames_grouped{area_i} = {param_grouped{area_i}.filename};
    
    laser_grouped{area_i} = [param_grouped{area_i}.laser_V];
%     laser_grouped{area_i}(laser_grouped{area_i}>5) = 5;
    nCellsByArea{area_i} = numel(filenames_grouped{area_i});
    
    grouped_trace{area_i} = cat(3,EPSC_grouped{area_i}.meanTrace);
    denominator = min(grouped_trace{area_i}(:,2,:));
    temptrace = grouped_trace{area_i}./abs(denominator);
    norm_traces{area_i} = squeeze(temptrace);
    
    collapse_minCurrent = [collapse_minCurrent; vertcat(EPSC_grouped{area_i}.minCurrent)];
%     collapse_minCurrent = [collapse_minCurrent; vertcat(EPSC_grouped{area_i}.tMinCurrentMean)];
    
    temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceRise);
    temp_cat(temp_cat>1) = temp_cat(temp_cat>1)/20000;
%     temp_cat(temp_cat>0.003) = NaN;
    pyr_rise_time{area_i} = temp_cat(:,2)*1000;
    IN_rise_time{area_i} = temp_cat(:,1)*1000;
    
    temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceLatency);
    pyr_latency{area_i} = temp_cat(:,2)./20;
    IN_latency{area_i} = temp_cat(:,1)./20;
    
    temp_cat = vertcat(EPSC_grouped{area_i}.jitter);

    if strcmp(cell_name,'SOM')
        IN_jitter{area_i} = temp_cat(:,1)*1000*.9;
    else
        IN_jitter{area_i} = temp_cat(:,1)*1000;
    end
    pyr_jitter{area_i} = temp_cat(:,2)*1000;
end
minCurrent = {collapse_minCurrent(:,1) collapse_minCurrent(:,2)};
timeVector = make_time(grouped_trace{area_i},20000,1);
temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);

traces_lat = cat(3,grouped_trace{1},grouped_trace{2});
traces_med = cat(3,grouped_trace{3},grouped_trace{4});


traces_pos = cat(3,grouped_trace{1},grouped_trace{3});
traces_ant = cat(3,grouped_trace{2},grouped_trace{4});

% calculate early and late EPSCs

for area_i = 1:4
    for file_i = 1:numel(EPSC_grouped{area_i})
        [temp_denom,temp_idx] = min(EPSC_grouped{area_i}(file_i).meanTrace(30:250,:,:));
        temp_idx = temp_idx + 30;
        if setMin
            temp_idx = EPSC_grouped{area_i}(file_i).earlyLatePoint;
        end
            temp_lat = EPSC_grouped{area_i}(file_i).meanTraceLatency;
            temp_trace_norm{area_i}(file_i,:,1) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(1):temp_lat(1)+200,1)),0,1);
            temp_trace_norm{area_i}(file_i,:,2) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(2):temp_lat(2)+200,2)),0,1);

            temp_trace{area_i}(file_i,:,:) = EPSC_grouped{area_i}(file_i).meanTrace;
        
        IN_early_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,1:temp_idx(1),1));
        IN_late_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,temp_idx(1):end,1));
        
        pyr_early_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,1:temp_idx(2),2));
        pyr_late_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,temp_idx(2):end,2));
        
        IN_early_EPSC_norm{area_i}(file_i) = (temp_trace_norm{area_i}(file_i,temp_idx(1),1));
        IN_late_EPSC_norm{area_i}(file_i) = 1-(temp_trace_norm{area_i}(file_i,temp_idx(1),1));
        
        pyr_early_EPSC_norm{area_i}(file_i) = (temp_trace_norm{area_i}(file_i,temp_idx(2),2));
        pyr_late_EPSC_norm{area_i}(file_i) = 1-(temp_trace_norm{area_i}(file_i,temp_idx(2),2));
        
        half_point_IN{area_i}(file_i) = find(temp_trace_norm{area_i}(file_i,:,1)>0.5,1,'first')./20;
        half_point_pyr{area_i}(file_i) = find(temp_trace_norm{area_i}(file_i,:,2)>0.5,1,'first')./20;
    end
end

nCellsByArea = cellfun(@numel,IN_early_EPSC,'un',0);

temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});

[p_IN_early_EPSC,~,stats_result_IN_early_EPSC] = anova1(cell2mat(IN_early_EPSC)',area_ids,'off');
% multcompare(stats_result_IN_early_EPSC)

[p_IN_late_EPSC,~,stats_result_IN_late_EPSC] = anova1(cell2mat(IN_late_EPSC)',area_ids,'off');
% multcompare(stats_result_IN_late_EPSC)

[p_pyr_early_EPSC,~,stats_result_pyr_early_EPSC] = anova1(cell2mat(pyr_early_EPSC)',area_ids,'off');
% multcompare(stats_result_pyr_early_EPSC)

[p_pyr_late_EPSC,~,stats_result_pyr_late_EPSC] = anova1(cell2mat(pyr_late_EPSC)',area_ids,'off');
% multcompare(stats_result_pyr_late_EPSC)

[p_IN_early_EPSC_norm,~,stats_result_IN_early_EPSC_norm] = anova1(cell2mat(IN_early_EPSC_norm)',area_ids,'off');
% multcompare(stats_result_IN_early_EPSC_norm)

[p_IN_late_EPSC_norm,~,stats_result_IN_late_EPSC_norm] = anova1(cell2mat(IN_late_EPSC_norm)',area_ids,'off');
% multcompare(stats_result_IN_late_EPSC_norm)

[p_pyr_early_EPSC_norm,~,stats_result_pyr_early_EPSC_norm] = anova1(cell2mat(pyr_early_EPSC_norm)',area_ids,'off');
% multcompare(stats_result_pyr_early_EPSC_norm)

[p_pyr_late_EPSC_norm,~,stats_result_pyr_late_EPSC_norm] = anova1(cell2mat(pyr_late_EPSC_norm)',area_ids,'off');
% multcompare(stats_result_pyr_late_EPSC_norm)
%%
fs = 10;
temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});
figure;
% onset latency
subplot(3,2,1);
plotSpread(pyr_latency,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,pyr_latency,1,'continuous',false,'cells_as_x',true);
fix_axes(gcf,fs,'Area','exc latency');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
axis square;
subplot(3,2,2);
plotSpread(IN_latency,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,IN_latency,1,'continuous',false);
fix_axes(gcf,fs,'Area','inh latency');
ylim([2 6]); 
axis square;

% jitter
% onset latency
subplot(3,2,3);
plotSpread(pyr_jitter,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,pyr_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area','exc jitter');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
ylim([0 2]);
axis square;

subplot(3,2,4);
plotSpread(IN_jitter,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,IN_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area','inh jitter');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
ylim([0 2]);
axis square;

subplot(3,2,5);
plotSpread(pyr_rise_time,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,pyr_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area','exc rise');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
ylim([0 2]);
axis square;

subplot(3,2,6);
plotSpread(IN_rise_time,'distributionMarkers','o','binWidth',1);
fast_errbar(1:4,IN_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area','inh rise');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
ylim([0 2]);
axis square;


%%
doZ = true;
if doZ
    IN_jitter_new = nanscore(cell2mat(IN_jitter'));
    IN_latency_new = zscore(cell2mat(IN_latency'));
    IN_amp_new = zscore(abs(collapse_minCurrent(:,1)));

    pyr_jitter_new = nanscore(cell2mat(pyr_jitter'));
    pyr_latency_new = zscore(cell2mat(pyr_latency'));
    pyr_amp_new = zscore(abs(collapse_minCurrent(:,2)));
    

    grouped_features = [{[pyr_jitter_new pyr_latency_new pyr_amp_new]};
        {[IN_jitter_new IN_latency_new IN_amp_new]}];
    
    cluster_ids = cellfun(@(x) kmeans(x(:,[1 2]),2),grouped_features,'un',0);
else
    IN_jitter_new = (cell2mat(IN_jitter'));
    IN_latency_new = (cell2mat(IN_latency'));
    IN_amp_new = abs(collapse_minCurrent(:,1));

    pyr_jitter_new = (cell2mat(pyr_jitter'));
    pyr_latency_new = (cell2mat(pyr_latency'));
    pyr_amp_new = abs(collapse_minCurrent(:,2));
    
grouped_features = [{[pyr_jitter_new pyr_latency_new pyr_amp_new]};
    {[IN_jitter_new IN_latency_new IN_amp_new]}];
    
    cluster_ids = cellfun(@(x) kmeans(x(:,[1 2]),2),grouped_features,'un',0);

end


cluster_ids_byArea = mat2cell(cell2mat(cluster_ids')',2,cell2mat(cellfun(@numel,pyr_jitter,'un',0)));

%%
fs = 10;
figure; hold on;
for cluster_i = 1:2
    subplot(2,2,cluster_i);
    plot_cluster = cellfun(@(x,y) x(logical(prod(y==cluster_i))),int_pyr_ratio_grouped,cluster_ids_byArea,'un',0);
    plotSpread(plot_cluster,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
    fast_errbar(1:4,plot_cluster,2,'continuous',false);
    axis square;
    xticks(1:4);
    xticklabels(area_names');
    xlim([0 5]);ylim([0 1.5]);
    fix_axes(gcf,fs,'Area','IN:Pyr');
    title(['matched pyr cluster',num2str(cluster_i),' IN cluster ',num2str(cluster_i)]); 
    
    subplot(2,2,cluster_i+2);
    plot_cluster = cellfun(@(x,y) x(logical((y(2,:)~=cluster_i).*(y(1,:)==cluster_i))),int_pyr_ratio_grouped,cluster_ids_byArea,'un',0);
    
    plotSpread(plot_cluster,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
    fast_errbar(1:4,plot_cluster,2,'continuous',false);
    axis square;
    xticks(1:4);ylim([0 1.5]);
    xticklabels(area_names');
    xlim([0 5]);
    fix_axes(gcf,fs,'Area','IN:Pyr');
    title(['unmatched pyr cluster',num2str(cluster_i),' IN cluster ',num2str(abs(cluster_i-3))]); 

end

%% PV and pyr by pyr grouping
fs = 10; ms = 4;

cluster1_colors = {[0 0 0],int_color};
cluster2_colors = cellfun(@(x) rescale(x,0.8,1),cluster1_colors,'un',0);

figure;

subplot(2,3,1);hold on;
cellfun(@(x,y,z) plot(x(y == 1,2),x(y == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,2),x(y == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Jitter'); axis square;
% ylim([0 2]);

subplot(2,3,2);
hold on;
cellfun(@(x,y,z) plot(x(y == 1,2),x(y == 1,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,2),x(y == 2,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Amplitude'); axis square;
% ylim([0 2]);

subplot(2,3,3);
hold on;
cellfun(@(x,y,z) plot(x(y == 1,3),x(y == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,3),x(y == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Amplitude','Jitter'); axis square;

subplot(2,3,4);hold on;
cellfun(@(x,z) plot(x(cluster_ids{1} == 1,2),x(cluster_ids{1} == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster1_colors','un',0);
cellfun(@(x,z) plot(x(cluster_ids{1} == 2,2),x(cluster_ids{1} == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Jitter'); axis square;
% ylim([0 2]);

subplot(2,3,5);
hold on;
cellfun(@(x,z) plot(x(cluster_ids{1} == 1,2),x(cluster_ids{1} == 1,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster1_colors','un',0);
cellfun(@(x,z) plot(x(cluster_ids{1} == 2,2),x(cluster_ids{1} == 2,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Amplitude'); axis square;
% ylim([0 2]);

subplot(2,3,6);
hold on;
cellfun(@(x,z) plot(x(cluster_ids{1} == 1,3),x(cluster_ids{1} == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster1_colors','un',0);
cellfun(@(x,z) plot(x(cluster_ids{1} == 2,3),x(cluster_ids{1} == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster2_colors','un',0);
fix_axes(gcf,fs,'Amplitude','Jitter'); axis square;
% ylim([0 2]);

%%
fs = 10; ms = 4;
figure;
% plot by cell types
subplot(2,3,1);
hold on;
plot((pyr_latency_new'),(pyr_jitter_new'),'ko','MarkerFaceColor','k','MarkerSize',ms);
plot((IN_latency_new'),(IN_jitter_new'),'o','Color',int_color,'MarkerFaceColor',int_color,'MarkerSize',ms);
fix_axes(gcf,fs,'Latency','Jitter'); axis square;

subplot(2,3,2);
hold on;hold on;
plot((pyr_latency_new'),(pyr_amp_new'),'ko','MarkerFaceColor','k','MarkerSize',ms);
plot((IN_latency_new'),(IN_amp_new'),'o','Color',int_color,'MarkerFaceColor',int_color,'MarkerSize',ms);
fix_axes(gcf,fs,'Latency','Amplitude'); axis square;
 
subplot(2,3,3);
hold on;
hold on;
plot((pyr_amp_new'),(pyr_jitter_new'),'ko','MarkerFaceColor','k','MarkerSize',ms);
plot((IN_amp_new'),(IN_jitter_new'),'o','Color',int_color,'MarkerFaceColor',int_color,'MarkerSize',ms);
fix_axes(gcf,fs,'Amplitude','Jitter'); axis square;


cluster1_colors = {[0 0 0],int_color};
cluster2_colors = cellfun(@(x) rescale(x,0.8,1),cluster1_colors,'un',0);

subplot(2,3,4);hold on;
cellfun(@(x,y,z) plot(x(y == 1,2),x(y == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,2),x(y == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Jitter'); axis square;
% ylim([0 2]);

subplot(2,3,5);
hold on;
cellfun(@(x,y,z) plot(x(y == 1,2),x(y == 1,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,2),x(y == 2,3),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Latency','Amplitude'); axis square;
% ylim([0 2]);

subplot(2,3,6);
hold on;
cellfun(@(x,y,z) plot(x(y == 1,3),x(y == 1,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster1_colors','un',0);
cellfun(@(x,y,z) plot(x(y == 2,3),x(y == 2,1),'o','Color',z,'MarkerFaceColor',z,'MarkerSize',ms),grouped_features,cluster_ids,cluster2_colors','un',0);
fix_axes(gcf,fs,'Amplitude','Jitter'); axis square;
% ylim([0 2]);

%%
figure;
subplot(1,2,1);
plot(cell2mat(pyr_latency'),cell2mat(pyr_jitter'),'o');
axis square;
ylim([0 2]);
xlim([2 5]);
fix_axes(gcf,20,'Latency (ms)','Jitter (ms)');
subplot(1,2,2);
plot(cell2mat(IN_latency'),cell2mat(IN_jitter'),'o');
axis square;
ylim([0 2]);
xlim([2 5]);
fix_axes(gcf,20,'Latency (ms)','Jitter (ms)');

%%

figure;
subplot(1,2,1);
temp=abs((collapse_minCurrent'));
plot(temp(2,:),cell2mat(pyr_jitter'),'o');
axis square;
ylim([0 2]);
% xlim([2 5]);
fix_axes(gcf,20,'pA','Jitter (ms)');
subplot(1,2,2);

plot(temp(1,:),cell2mat(IN_jitter'),'o');
axis square;
ylim([0 2]);
% xlim([2 5]);
fix_axes(gcf,20,'pA','Jitter (ms)');

%%
all_hvas = {'LM','AL','PM','AM'};
timeVector = 1000*make_time(temp_trace{1},20000,2);

figure;
subplot(2,3,1)
plotSpread(IN_early_EPSC,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_early_EPSC,2,'continuous',false);
axis square; 
title([cell_name,' early EPSC']);
xticks(1:4);xticklabels(all_hvas); 

subplot(2,3,2)
plotSpread(IN_late_EPSC,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_late_EPSC,2,'continuous',false);
axis square;
title([cell_name,' late EPSC']);
xticks(1:4);xticklabels(all_hvas); 

subplot(2,3,3); hold on;
for area_i = 1:4
    plot(timeVector,mean(temp_trace{area_i}(:,:,1)),'Color',cmap(area_i,:));
end
axis square;axis tight;
title([cell_name,' avg traces']);

subplot(2,3,4)
plotSpread(pyr_early_EPSC,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_early_EPSC,2,'continuous',false);
axis square; 
title('pyr early EPSC');
xticks(1:4);xticklabels(all_hvas); 
subplot(2,3,5)
plotSpread(pyr_late_EPSC,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_late_EPSC,2,'continuous',false);
axis square;
title('pyr late EPSC');
xticks(1:4);xticklabels(all_hvas); 

subplot(2,3,6); hold on;
for area_i = 1:4
    plot(timeVector,mean(temp_trace{area_i}(:,:,2)),'Color',cmap(area_i,:));
end
axis square;axis tight;
title('pyr avg traces');

%%

all_hvas = {'LM','AL','PM','AM'};
fs = 10;
figure;
subplot(1,4,1)
plotSpread(half_point_IN,'distributionMarkers','o','binWidth',1,'distributionColors',cmap);
fast_errbar(1:4,half_point_IN,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','half point (s)');
title([cell_name,' half point']);
xticks(1:4);xticklabels(all_hvas); 
ylim([0 10]);

subplot(1,4,2)
plotSpread(half_point_pyr,'distributionMarkers','o','binWidth',1,'distributionColors',cmap);
fast_errbar(1:4,half_point_pyr,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','half point (s)');
title(['pyr half point']);
xticks(1:4);xticklabels(all_hvas);
ylim([0 10]);

figure;
[p_half_IN,~,stats_result_half_IN] = anova1(cell2mat(half_point_IN)',area_ids,'off');
multcompare(stats_result_half_IN)

figure;
[p_half_pyr,~,stats_result_half_pyr] = anova1(cell2mat(half_point_pyr)',area_ids,'off');
multcompare(stats_result_half_pyr)

%%
all_hvas = {'LM','AL','PM','AM'};
timeVector = 1000*make_time(temp_trace_norm{1},20000,2);
fs = 10;
figure;
subplot(2,3,1)
plotSpread(IN_early_EPSC_norm,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_early_EPSC_norm,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','norm charge');
title([cell_name,' early EPSC - norm.']);
xticks(1:4);xticklabels(all_hvas); ylim([0 1]);

subplot(2,3,2)
plotSpread(IN_late_EPSC_norm,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_late_EPSC_norm,2,'continuous',false);
axis square;fix_axes(gcf,fs,'Area','norm charge');
title([cell_name,' late EPSC - norm.']);ylim([0 1]);
xticks(1:4);xticklabels(all_hvas); 

subplot(2,3,3); hold on;
for area_i = 1:4
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,1)),'Color',cmap(area_i,:));
end
axis square;
title([cell_name,' mean trace - norm.']);
axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');

subplot(2,3,4)
plotSpread(pyr_early_EPSC_norm,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_early_EPSC_norm,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','norm charge');
title('pyr early EPSC - norm.');
xticks(1:4);xticklabels(all_hvas); ylim([0 1]);

subplot(2,3,5)
plotSpread(pyr_late_EPSC_norm,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_late_EPSC_norm,2,'continuous',false);
axis square;fix_axes(gcf,fs,'Area','norm charge');
title('pyr late EPSC - norm.');
xticks(1:4);xticklabels(all_hvas); ylim([0 1]);

subplot(2,3,6); hold on;
for area_i = 1:4
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,2)),'Color',cmap(area_i,:));
end
axis square;axis tight;
title('pyr avg traces - norm.');
fix_axes(gcf,fs,'time (ms)','cumulative charge');

%%
figure;
timeVector = make_time(temp_trace_norm{1},2,2);
fs = 10;
subplot(1,4,1); hold on;
for area_i = 1:4
%     plot(timeVector,temp_trace_norm{area_i}(:,:,2),'Color',rescale(cmap(area_i,:),.8,1));
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,2)),'Color',cmap(area_i,:));
end
axis square;axis tight;
title('pyr avg traces - norm.');
fix_axes(gcf,fs,'time (ms)','cumulative charge');

subplot(1,4,2); hold on;
for area_i = 1:4
%     plot(timeVector,temp_trace_norm{area_i}(:,:,1),'Color',rescale(cmap(area_i,:),.8,1));
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,1)),'Color',cmap(area_i,:));
end
axis square;
title([cell_name,' mean trace - norm.']);
axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');

%% latency, early late, rise time, jitter
fs = 9;
diff_rise_time = cellfun(@(x,y) x-y,IN_rise_time,pyr_rise_time,'un',0);

temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});

temp_earlyLate = cellfun(@(x) vertcat(x.earlyLatePoint),EPSC_grouped,'un',0);
IN_earlyLate = cellfun(@(x) x(:,1)./20,temp_earlyLate,'un',0);
pyr_earlyLate = cellfun(@(x) x(:,2)./20,temp_earlyLate,'un',0);

[p_half_INearlylate,~,stats_result_INearlylate] = anova1(cell2mat(pyr_earlyLate)',area_ids,'off');
multcompare(stats_result_INearlylate)

% quick stats
[p_pyr_rise_time,~,stats_result_pyr_rise_time] = anova1(cell2mat(pyr_rise_time'),area_ids,'off');
multcompare(stats_result_pyr_rise_time)

[p_IN_rise_time,~,stats_result_IN_rise_time] = anova1(cell2mat(IN_rise_time'),area_ids,'off');
multcompare(stats_result_IN_rise_time)


clear ax
figure;
subplot(1,4,1);
plotSpread(IN_latency,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_latency,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' latency']); axis square; xticks(1:4);xticklabels(area_names);ylim([0 8]);yticks(0:2:10);
xlim([0 5]);

subplot(1,4,2);
plotSpread(IN_earlyLate,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap)
fast_errbar(1:4,IN_earlyLate,1,'continuous',false);
axis square; xticks(1:4); xticklabels(area_names); ylim([0 8]);yticks(0:2:10);
fix_axes(gcf,fs,'Area',[cell_name, ' early point']);

subplot(1,4,3);
plotSpread(IN_rise_time,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' 10-90 rise time']); axis square; xticks(1:4);xticklabels(area_names);
ylim([0 2]);

subplot(1,4,4);
plotSpread(IN_jitter,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' IN jitter']); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:0.2:1.2); 
yticks([0:0.5:1.5]); ylim([0 1.5]);
xlim([0 5]);
%%
subplot(1,4,1);
plotSpread(pyr_latency,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_latency,1,'continuous',false);
fix_axes(gcf,fs,'Area','pyr 10-90 latency'); axis square; xticks(1:4);xticklabels(area_names);yticks(0:2:10);ylim([0 8]);
% 
subplot(1,4,2);
plotSpread(pyr_earlyLate,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap)
fast_errbar(1:4,pyr_earlyLate,1,'continuous',false); 
axis square;fix_axes(gcf,fs,'Area',['pyr early late point (ms)'])
xticks(1:4); xticklabels(area_names); yticks(0:2:10);ylim([0 8]); 

subplot(1,4,3);
plotSpread(pyr_rise_time,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area','pyr 10-90 rise time'); axis square; xticks(1:4);xticklabels(area_names);
ylim([0 2]);

subplot(1,4,4);
plotSpread(pyr_jitter,'spreadWidth',0.75,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,pyr_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area','pyr jitter'); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:0.2:0.6);ylim([0 0.6]);


%% single cell early late

figure; subplot(1,3,1); hold on;
cellfun(@(x,y,z) plot(x,y,'o','Color',z),pyr_latency,IN_latency,cmap_cell','un',0);
fix_axes(gcf,15,'pyr latency (ms)',[cell_name,' latency (ms)']); axis square; matchxy('min');
subplot(1,3,2); hold on;
cellfun(@(x,y,z) plot(x,y,'o','Color',z),pyr_jitter,IN_jitter,cmap_cell','un',0);
fix_axes(gcf,15,'pyr jitter (ms)',[cell_name,' jitter (ms)']); axis square; matchxy('min');
subplot(1,3,3); hold on;
cellfun(@(x,y,z) plot(x,y,'o','Color',z),pyr_earlyLate,IN_earlyLate,cmap_cell,'un',0);
fix_axes(gcf,15,'pyr earlyLate (ms)',[cell_name,' earlyLate (ms)']); axis square; matchxy('min');

%% latency, early late, rise time, jitter as difference
fs = 10;

temp_earlyLate = cellfun(@(x) vertcat(x.earlyLatePoint),EPSC_grouped,'un',0);
IN_earlyLate = cellfun(@(x) x(:,1)./20,temp_earlyLate,'un',0);
pyr_earlyLate = cellfun(@(x) x(:,2)./20,temp_earlyLate,'un',0);

diff_latency = cellfun(@(x,y) x-y,IN_latency,pyr_latency,'un',0);
diff_earlyLate = cellfun(@(x,y) x-y,IN_earlyLate,pyr_earlyLate,'un',0);
diff_rise_time = cellfun(@(x,y) x-y,IN_rise_time,pyr_rise_time,'un',0);
diff_jitter = cellfun(@(x,y) x-y,IN_jitter,pyr_jitter,'un',0);

% quick stats
[p_pyr_rise_time,~,stats_result_pyr_rise_time] = anova1(cell2mat(diff_rise_time'),area_ids,'off');
multcompare(stats_result_pyr_rise_time)

[p_IN_rise_time,~,stats_result_IN_rise_time] = anova1(cell2mat(IN_rise_time'),area_ids,'off');
multcompare(stats_result_IN_rise_time)


clear ax
figure;
subplot(2,4,1);
plotSpread(IN_latency,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_latency,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' 10-90 latency']); axis square; xticks(1:4);xticklabels(area_names);ylim([0 10]);yticks(0:2:10);

subplot(2,4,2);
plotSpread(IN_earlyLate,'distributionMarkers','o','distributionColors',cmap)
fast_errbar(1:4,IN_earlyLate,1,'continuous',false);
axis square; xticks(1:4); xticklabels(area_names); ylim([0 10]);yticks(0:2:10);
fix_axes(gcf,fs,'Area',[cell_name, ' early late point (ms)']);

subplot(2,4,3);
plotSpread(IN_rise_time,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' 10-90 rise time']); axis square; xticks(1:4);xticklabels(area_names);
ylim([0 3]);
subplot(2,4,4);
plotSpread(IN_jitter,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,IN_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area',[cell_name,' IN jitter']); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:0.2:0.6);ylim([0 0.6]);

subplot(2,4,5);
plotSpread(diff_latency,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,diff_latency,1,'continuous',false);
fix_axes(gcf,fs,'Area','diff with pyr 10-90 latency'); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:2:10);ylim([0 10]);

subplot(2,4,6);
plotSpread(diff_earlyLate,'distributionMarkers','o','distributionColors',cmap)
fast_errbar(1:4,diff_earlyLate,1,'continuous',false); 
axis square;fix_axes(gcf,fs,'Area','diff with pyr early late point (ms)')
xticks(1:4); xticklabels(area_names); 
% yticks(0:2:10);ylim([0 10]); 

subplot(2,4,7);
plotSpread(diff_rise_time,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,diff_rise_time,1,'continuous',false);
fix_axes(gcf,fs,'Area','diff with pyr 10-90 rise time'); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:4);

subplot(2,4,8);
plotSpread(diff_jitter,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,diff_jitter,1,'continuous',false);
fix_axes(gcf,fs,'Area','diff with pyr jitter'); axis square; xticks(1:4);xticklabels(area_names);
% yticks(0:0.2:0.6);ylim([0 0.6]);


%% plot grand average in 4 areas int EPSC vs Pyr EPSC
fs = 10;
clear ax
doNorm = false;
normLM = true;
LM_denom = -min(mean(grouped_trace{1},3));
for area_i = 1:4
    if doNorm
        
        figure(1); hold on;
        ax(area_i) = subplot(1,4,area_i); hold on;
        fast_errbar(timeVector,norm_traces{area_i}(:,1,:),3,'shaded',true,'Color',int_color);
        fast_errbar(timeVector,norm_traces{area_i}(:,2,:),3,'shaded',true,'Color',[0 0 0]);
        title(area_names{area_i});
        fix_axes(gcf,10,'time (ms)','norm pA');
        axis square
    elseif normLM
        figure(1); hold on;
        
        ax(area_i) = subplot(1,4,area_i); hold on;
        fast_errbar(timeVector,grouped_trace{area_i}(:,1,:)./LM_denom(1),3,'shaded',true,'Color',int_color);
        fast_errbar(timeVector,grouped_trace{area_i}(:,2,:)./LM_denom(2),3,'shaded',true,'Color',[0 0 0]);
        title(area_names{area_i});
        fix_axes(gcf,10,'time (ms)','pA');
        axis square
    else
        figure(1); hold on;
        ax(area_i) = subplot(1,4,area_i); hold on;
        fast_errbar(timeVector,grouped_trace{area_i}(:,1,:),3,'shaded',true,'Color',int_color);
        fast_errbar(timeVector,grouped_trace{area_i}(:,2,:),3,'shaded',true,'Color',[0 0 0]);
        title(area_names{area_i});
        fix_axes(gcf,10,'time (ms)','pA');
        axis square
    end
    
    
end
linkaxes(ax);

figure(2); hold on;
subplot(1,4,1);
plotSpread(int_pyr_ratio_grouped,'binWidth',1,'distributionMarkers','o','distributionColors',cmap);
fast_errbar(1:4,int_pyr_ratio_grouped,2,'continuous',false,'stats',true);
axis square;
xticks(1:4);
xticklabels(area_names');
xlim([0 5]);
fix_axes(gcf,fs,'Area','IN:Pyr');
% axis square
%%
[p,~,stats_result] = anova1((cell2mat(int_pyr_ratio_grouped)'),area_ids);
% [p,~,stats_result] = kruskalwallis((cell2mat(int_pyr_ratio_grouped)'),area_ids);
multcompare(stats_result)

%%

early_sumMed = squeeze(sum(traces_med(1:earlyLate_bound,:,:)));
early_sumLat = squeeze(sum(traces_lat(1:earlyLate_bound,:,:)));

%% plot grand average in medial vs lateral areas for int and pyr neurons
clear ax
figure; 
doNorm = true;
medLat = true;

if medLat
    traces_1 = traces_lat;
    traces_2 = traces_med;
    titles = [{'lateral'},{'medial'}];
else
    traces_1 = traces_ant;
    traces_2 = traces_pos;
    titles = [{'anterior'},{'posterior'}];
end

if doNorm
    ax(1) = subplot(1,2,1);hold on;
    fast_errbar([],traces_1(:,1,:)./abs(min(traces_1(30:120,2,:))),3,'shaded',true,'Color',int_color)
    fast_errbar([],traces_1(:,2,:)./abs(min(traces_1(30:120,2,:))),3,'shaded',true,'Color',[0 0 0])
    fix_axes(gcf,20,'time (ms)','pA');
    title(titles{1});
    axis square
    ax(2) = subplot(1,2,2);hold on;
    fast_errbar([],traces_2(:,1,:)./abs(min(traces_2(30:120,2,:))),3,'shaded',true,'Color',int_color)
    fast_errbar([],traces_2(:,2,:)./abs(min(traces_2(30:120,2,:))),3,'shaded',true,'Color',[0 0 0])
    title(titles{2});
    axis square
    linkaxes(ax);
    fix_axes(gcf,20,'time (ms)','pA');
else
    ax(1) = subplot(1,2,1);hold on;
    fast_errbar([],traces_1(:,1,:),3,'shaded',true,'Color',int_color)
    fast_errbar([],traces_1(:,2,:),3,'shaded',true,'Color',[0 0 0])
    fix_axes(gcf,20,'time (ms)','pA');
    title(titles{1});
    axis square
    ax(2) = subplot(1,2,2);hold on;
    fast_errbar([],traces_2(:,1,:),3,'shaded',true,'Color',int_color)
    fast_errbar([],traces_2(:,2,:),3,'shaded',true,'Color',[0 0 0])
    title(titles{2});
    axis square
    linkaxes(ax);
    fix_axes(gcf,20,'time (ms)','pA');
end
% plot ratio
figure;
hold on;
plot(1,cell2mat(int_pyr_ratio_grouped(1:2)),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',20);
plot(1,mean(cell2mat(int_pyr_ratio_grouped(1:2))),'ko','LineWidth',5,'MarkerSize',20);
plot(2,cell2mat(int_pyr_ratio_grouped(3:4)),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',20);
plot(2,mean(cell2mat(int_pyr_ratio_grouped(3:4))),'ko','LineWidth',5,'MarkerSize',20);
xticks(1:2);xticklabels({'Lateral','Medial'});xlim([0 3]);
fix_axes(gcf,40,'Area','IN:Pyr');
%%
% calculate statistics
lat_med_id = ones(1,sum(cell2mat(nCellsByArea)));
lat_med_id(logical(sum(area_ids == 'L',2))) = 0;

temp = cell2mat(int_pyr_ratio_grouped);
ranksum(temp(logical(sum(area_ids == 'L',2))),temp(~logical(sum(area_ids == 'L',2))))

%% compare med/lat int EPSCS
doNorm= true;
normSelf = true;
medLat = true;

if medLat
    traces_1 = traces_lat;
    traces_2 = traces_med;
    titles = [{'lateral'},{'medial'}];
else
    traces_1 = traces_ant;
    traces_2 = traces_pos;
    titles = [{'anterior'},{'posterior'}];
end
figure;hold on;
for cell_type = 1:2
    subplot(1,2,cell_type); hold on;
    color = [int_color;0.5 0.5 0.5];
    if doNorm
        if normSelf
            fast_errbar([],squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,cell_type,:)))),2,'shaded',true,'color',(color(cell_type,:)-0.3));
            fast_errbar([],squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,cell_type,:)))),2,'shaded',true,'color',color(cell_type,:));
        else
            fast_errbar([],squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,2,:)))),2,'shaded',true,'color',(color(cell_type,:)-0.3))
            fast_errbar([],squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,2,:)))),2,'shaded',true,'color',(color(cell_type,:)))
        end
    else
        fast_errbar([],squeeze(traces_1(:,cell_type,:)),2,'shaded',true,'color',(color(cell_type,:)-0.3))
        fast_errbar(squeeze(traces_2(:,cell_type,:)),2,'shaded',true,'color',(color(cell_type,:)))
    end
    legend([titles(2);{''};titles(1);{''}]);
    axis tight; axis square;
end

%%
timeVector = make_time(grouped_trace{area_i},20000,1);
area_i = 1;
cell_i = 1;

while cell_i <= nCellsByArea{area_i}
figure(100);clf;
subplot(2,2,1);hold on;

for trial_i = 1:5
plot(timeVector,smooth(squeeze(EPSC_grouped{area_i}(cell_i).trace(:,2,trial_i))),'Color',[0.8 0.8 0.8]);
plot(timeVector,smooth(squeeze(EPSC_grouped{area_i}(cell_i).trace(:,1,trial_i))),'Color',[1 0.8 0.8]);    
end

plot(timeVector,smooth(mean(squeeze(EPSC_grouped{area_i}(cell_i).trace(:,2,:)),2)),'k');
plot(timeVector,smooth(mean(squeeze(EPSC_grouped{area_i}(cell_i).trace(:,1,:)),2)),'r');
vline(EPSC_grouped{area_i}(cell_i).earlyLatePoint(1)/20000,'r--');
vline(EPSC_grouped{area_i}(cell_i).earlyLatePoint(2)/20000,'k--');
fix_axes(gcf,10,'time (s)','pA');
title(['IN:Pyr:',num2str(int_pyr_ratio_grouped{area_i}(cell_i))]);
axis square;
subplot(2,2,2);hold on;
plot(timeVector(30:120),squeeze(EPSC_grouped{area_i}(cell_i).trace(30:120,2,:)),'Color',[0.8 0.8 0.8]);
plot(timeVector(30:120),mean(squeeze(EPSC_grouped{area_i}(cell_i).trace(30:120,2,:)),2),'k');
plot(timeVector(30:120),squeeze(EPSC_grouped{area_i}(cell_i).trace(30:120,1,:)),'Color',[1 0.8 0.8]);
plot(timeVector(30:120),mean(squeeze(EPSC_grouped{area_i}(cell_i).trace(30:120,1,:)),2),'r');
axis tight;ylim([-50 10]);
fix_axes(gcf,10,'time(s)','pA');
axis square;

disp(EPSC_grouped{area_i}(cell_i).meanTraceRise);
% 
% subplot(2,2,3); hold on;
% plot(timeVector(1:end-1),(smooth(diff(smooth(EPSC_grouped{area_i}(cell_i).meanTrace(:,2))),10)),'k');
% plot(timeVector(1:end-1),(smooth(diff(smooth(EPSC_grouped{area_i}(cell_i).meanTrace(:,1))),10)),'r');
% axis square;
% 
% subplot(2,2,4);hold on;
% plot(rescale(cumsum(abs(EPSC_grouped{area_i}(cell_i).meanTrace(:,2))),0,1),'k');
% plot(rescale(cumsum(abs(EPSC_grouped{area_i}(cell_i).meanTrace(:,1))),0,1),'r');

wait_flag = input('next cell = enter; prev cell = 0');

if isempty(wait_flag)
    cell_i = cell_i + 1;
elseif wait_flag == 0
    cell_i = cell_i -1;
end
end
%%
fs = 15;
figure; 
subplot(1,2,1); hold on;
for area_i = 1:4
    tempY = abs(vertcat(EPSC_grouped{area_i}.minCurrent));
%     tempY = abs([IN_early_EPSC{area_i}' pyr_early_EPSC{area_i}']);
    plot(tempY(:,2),tempY(:,1),'o','color',cmap(area_i,:),'MarkerFaceColor',cmap(area_i,:),'MarkerSize',4);
    result = fitlm(tempY(:,2),tempY(:,1),'linear','Intercept',false);
    slope(area_i) = result.Coefficients.Estimate;
    slope_CI(area_i,:) = result.coefCI;
    
end
xvals = linspace(0,1000,100);
for area_i = 1:4
    plot(xvals,xvals*slope(area_i),'-','Color',cmap(area_i,:));
%     plot(xvals,xvals*slope_CI(area_i,1),'--','Color',cmap(area_i,:));
%     plot(xvals,xvals*slope_CI(area_i,2),'--','Color',cmap(area_i,:));
end
fix_axes(gcf,fs,'Pyr EPSC (pA)',[cell_name, ' EPSC (pA)']); axis square;
xlim([0 1000]); ylim([0 4000]);



%%
subplot(1,3,2); hold on;
for area_i = 1:4
    tempY = abs(vertcat(EPSC_grouped{area_i}.minCurrent));
temp_group = [temp_group; tempY(:,2), IN_late_EPSC_norm{area_i}'];    
plot(tempY(:,2),IN_late_EPSC_norm{area_i},'o','color',cmap(area_i,:),'MarkerFaceColor',cmap(area_i,:),'MarkerSize',4);
tempY = abs(pyr_early_EPSC{area_i}');
    temp_group = [temp_group; tempY, IN_late_EPSC_norm{area_i}'];
    plot(tempY,IN_late_EPSC_norm{area_i},'o','color',cmap(area_i,:),'MarkerFaceColor',cmap(area_i,:),'MarkerSize',4);
end
fix_axes(gcf,fs,'Pyr EPSC (pA)',[cell_name, ' norm late charge']); axis square;
ylim([0 1]);

subplot(1,3,3);
plot(temp_group(:,1),temp_group(:,2),'.'); lsline
fix_axes(gcf,fs,'Pyr EPSC (pA)',[cell_name, ' norm late charge']); axis square;
ylim([0 1]);

%% segment out top percentile of late charges

tenthSize = floor(sum(cell2mat(nCellsByArea))/10);
[~,top10] = maxk(cell2mat(IN_late_EPSC_norm),tenthSize);
[~,bot10] = mink(cell2mat(IN_late_EPSC_norm),tenthSize);

IN_trace = cell2mat(cellfun(@(x) squeeze(x(:,1,:)),grouped_trace,'un',0));

timeVector = make_time(temp_trace{1},20000,2);
figure;
subplot(1,2,1);  hold on;
fast_errbar(timeVector*1000,-IN_trace(:,top10)./min(IN_trace(:,top10)),2,'shaded',true,'Color',int_color);
fast_errbar(timeVector*1000,-IN_trace(:,bot10)./min(IN_trace(:,bot10)),2,'shaded',true,'Color',int_color);
fix_axes(gcf,15,'time (ms)','norm. pA'); axis square;

[~,top10] = maxk(cell2mat(pyr_late_EPSC_norm),tenthSize);
[~,bot10] = mink(cell2mat(pyr_late_EPSC_norm),tenthSize);
pyr_trace = cell2mat(cellfun(@(x) squeeze(x(:,2,:)),grouped_trace,'un',0));

subplot(1,2,2);  hold on;
fast_errbar(timeVector*1000,-pyr_trace(:,top10)./min(pyr_trace(:,top10)),2,'shaded',true);
fast_errbar(timeVector*1000,-pyr_trace(:,bot10)./min(pyr_trace(:,bot10)),2,'shaded',true);
fix_axes(gcf,15,'time (ms)','norm. pA'); axis square;