
icell = best_cell_list(5);
trace_no_ad_avg = []; trace_cond_avg_750 = []; trace_cond_avg_250 = [];

for idelta = 1 : ndelta
    trace_no_ad_avg(idelta, :) = nanmean(trace_no_ad_merge_dfof{icell, idelta}, 1);
    trace_cond_avg_750(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 2}, 1);
    trace_cond_avg_250(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 1}, 1);
end

x = [length(trace_no_ad_avg), length(trace_cond_avg_750), length(trace_cond_avg_250)];
xmax = max(x);
y = [trace_no_ad_avg, trace_cond_avg_750, trace_cond_avg_250];
ymin = min(y(:));
ymax = max(y(:));

figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle_LL(num2str(icell))
for col = 1 : 3
    for idelta = 1 : ndelta
        subplot(ndelta, 3, col + 3*(idelta-1));
        if col == 1, plot(trace_no_ad_avg(idelta, :))
            if idelta == 1, title('no adapter'), end
%             if idelta == 8, text(0,1,['[',num2str(round(ymin,3)),', ', num2str(round(ymax,3)),']'],'Units','normalized'), end
        elseif col == 2, plot(trace_cond_avg_750(idelta, :))
             if idelta == 1, title('isi 750'), end
        elseif col == 3, plot(trace_cond_avg_250(idelta, :))
             if idelta == 1, title('isi 250'), end
        end
%         axis off
        xlim([0, xmax])
        ylim([ymin, ymax])
        xticks(0 : 30 : xmax)
        yticks(round(ymin*10)/10 : 0.1 : round(ymax*10)/10)
    end
end
saveas(gcf, ['dfof trace ', num2str(icell), '.jpg'])
close