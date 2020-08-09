
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

for col = 1 : 3
    figure('units','normalized','outerposition',[0 0 1/3 1]);
    for idelta = 1 : ndelta
        subplot(ndelta,1,idelta);
        if col == 1, plot(trace_no_ad_avg(idelta, :))
            if idelta == 1, title('no adapter'), end
            if idelta == 8, text(0,1,['[',num2str(round(ymin,3)),', ', num2str(round(ymax,3)),']'],'Units','normalized'), end
        elseif col == 2, plot(trace_cond_avg_750(idelta, :))
             if idelta == 1, title('isi 750'), end
        elseif col == 3, plot(trace_cond_avg_250(idelta, :))
             if idelta == 1, title('isi 250'), end
        end
        axis off
        xlim([0, xmax])
        ylim([ymin, ymax])   
    end
    savefig(['test', num2str(col)])
    close
end

%%
h1 = openfig('test1.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('test2.fig','reuse');
ax2 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h3 = figure; %create new figure
s1 = subplot(2,1,1); %create and get handle to the subplot axes
s2 = subplot(2,1,2);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
