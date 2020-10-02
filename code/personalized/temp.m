%% adp_ratio_a0t0_high_ad wrong

% wth = {}; wth2 = {};
for iset = 1:nset
    tt = squeeze(adp(iset).adp_ratio_a0t0_high_ad);
    wth2{iset} = tt(:,1)./tt(:,2)
    
    tt = squeeze(adp(iset).adp_ratio_a0t0);
    wth{iset} = tt(:,1)./tt(:,2)
end

%% mirror dfof_equiv_ad_targ wrt to resp_ad

tt = dfof_equiv_ad_targ - dfof_equiv_ad
sum(tt,2)

tt = dfof_equiv_ad_targ ./ dfof_equiv_ad - 1
sum(tt,2)

%% check from adp_across (mat generator)

dfof_equiv_ad_targ = squeeze(dfof_avg_merge(:, 8, 2:3)); %750-250
dfof_equiv_ad(ii, igap) = dfof_avg_ad(icell, idelta);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% trace grand average for area or for mouse

trace = struct;
for iset = 1:nset
    trace(iset).ncell = size(set{iset, 5}.trace_targ0_750, 1) - nlow(iset);
    trace(iset).ncell_neg = [sum(adp(iset).adp_neg{1}), sum(adp(iset).adp_neg{2})]; % why complimentary mirror?
    trace(iset).ncell_pos = [sum(adp(iset).adp_pos{1}), sum(adp(iset).adp_pos{2})];
    
    trace(iset).trace_avg{1} = nanmean(set{iset, 5}.trace_targ0_750(adp(iset).adp_pos{1},:),1);
    trace(iset).trace_avg{2} = nanmean(set{iset, 5}.trace_targ0_250(adp(iset).adp_pos{2},:),1);

%     trace(iset).trace_std{1,1} = nanstd(set{iset, 5}.trace_targ0_750(~low_ad_resp{iset},:),1);
%     trace(iset).trace_std{2,1} = nanstd(set{iset, 5}.trace_targ0_250(~low_ad_resp{iset},:),1);
end

%% trace for area

by_area_id = {[1,4,7], [2,5,8], [3,6]}; narea = length(by_area_id);

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
legend('V1', 'LM', 'LI', 'Location','northeast'); legend boxoff
% for iarea = 1 : narea
%     x = 1:trace_len; x2 = [x, fliplr(x)];
%     curve1 = trace_area_avg(iarea, 1:trace_len) + trace_area_std(iarea, 1:trace_len); curve2 = trace_area_avg(iarea, 1:trace_len) - trace_area_std(iarea, 1:trace_len);
%     inBetween = [curve1, fliplr(curve2)];
%     h = fill(x2, inBetween, color_list{iarea}, 'edgecolor','none'); 
%     h.FaceAlpha = 0.3;
% end

end
% saveas(gcf, ['trace across area isi ', num2str(igap)], 'jpg'); close

%%
for iarea = 1:3
    resp_ad(iarea) = max(trace_area_avg(iarea, 1:15)) - min(trace_area_avg(iarea, 1:15));
    resp_ad_targ(iarea) = max(trace_area_avg(iarea, 16:60)) - min(trace_area_avg(iarea, 16:40));
    % bug here: should not use min. should be frame-confined
    adp_area_trace(iarea) = resp_ad_targ(iarea)/resp_ad(iarea) - 1;
end
adp_area_trace % possible bug: too rough, use avg window if necessary


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
% saveas(gcf, ['trace across area-mouse isi ', num2str(igap)], 'jpg'); close

