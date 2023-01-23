size(trace_sem)

trace_avg_all = mean(mean(trace_sem, 1), 2);
trace_avg_all = squeeze(trace_avg_all);
plot(trace_avg_all)