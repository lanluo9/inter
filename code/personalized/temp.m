% tt = mean(mean(trace_avg, 2),1)

% tt = cellfun(@mean, trace_avg);
% t = squeeze(nanmean(tt, 2));

% tt = squeeze(sum(tc_align_ad,3));

% tt = isnan(npSub_tc);
% imagesc(tt)

tt = isnan(dfof_align_ad);
imshow3D(tt);
set(gcf, 'Position', get(0, 'Screensize'));

%%

t = sum(nansum(dfof_align_ad,3),2);
t = isnan(t);
find(t)

%%

tt = squeeze(dfof_align_ad(1103,:,:));
imagesc(tt)
%%
tt = squeeze(tc_align_ad(1103,:,:));
imagesc(tt)

%%

t = squeeze(npSub_tc(:,1103));

%%
ttt = sum(npSub_tc,1);
void_id = (ttt==0);
tc_corrected = npSub_tc(:, ~void_id);

