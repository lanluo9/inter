
histogram(ori_pref_cells, 8)
max(ori_pref_cells)

%%% bin 0-deg preferring cells by ad resp
clear idx
[count, idx] = histc(dfof_avg_ad(pref_0_cell_real),-0.01:0.05:ceil(max(dfof_avg_ad(pref_0_cell_real))*10)/10); % relax lower bound: one cell dfof = -0.0049
count = count(1:end-1);
count_nonzero_id = find(count~=0);
% histogram(dfof_avg_ad(pref_0_cell_real))
resp_bin_ad = accumarray(idx(:), dfof_avg_ad(pref_0_cell_real),[],@mean)
resp_bin_tg = [];
resp_bin_std = []; 
resp_bin_ste = zeros(length(count),2);
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell_real,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell_real,2),[],@mean);
resp_bin_std(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell_real,1),[],@std);
resp_bin_std(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell_real,2),[],@std);
resp_bin_ste(count~=0, :) = resp_bin_std(count~=0, :) ./ sqrt(count(count~=0)) % use std or ste?

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b.'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r.')
errorbar(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, resp_bin_std(:,1), 'b', 'LineStyle','none');
errorbar(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, resp_bin_std(:,2), 'r', 'LineStyle','none');
for itext = 1 : length(count_nonzero_id)
    text(resp_bin_ad(count_nonzero_id(itext)), ...
        resp_bin_tg(count_nonzero_id(itext),2)./resp_bin_ad(count_nonzero_id(itext)) + resp_bin_std(count_nonzero_id(itext),2) + 0.02, ...
        ['n=', num2str(count(count_nonzero_id(itext)))], 'HorizontalAlignment', 'center')
end
ylim([0,1])
legend('isi 250', 'isi 750')
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['Fig1E cells prefer 0 deg binned by ad-resp.jpg'])
% close
