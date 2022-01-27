% % Sample Data
% C = repmat({num2cell(1:11,1)},100,1);
% % Unpacked
% C2 = vertcat(C{:})
% 
% t2 = vertcat(dfof_ad_trial{:});
% histogram(t2);
% t3 = (t2);
% mean(t3)
% median(t3)

% tt = dfof_ad_trial{1,:};


median(vertcat(dfof_ad_trial{:})), median(vertcat(dfof_tg_trial{:}))
% dfof_ad
median(dfof_ad(:)), median(dfof_tg(:))

i = 15
median(median( (dfof_tg(:,i) - dfof_ad(:,i)) / (dfof_ad(:,i) + 1e-10) ))