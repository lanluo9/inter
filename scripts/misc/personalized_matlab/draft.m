tmp = load(bootstrap_file);
well_fit_bool = tmp.well_fit_cell;
sum(well_fit_bool) / length(well_fit_bool)

ori_pref_noad = ori_pref(:, 1); % first col is no-adapter ori pref
ori_pref_noad_boot_avg = mean(tmp.ori_pref_runs, 2); % avg across boots
ori_pref_noad_boot_med = median(tmp.ori_pref_runs, 2);

ori_pref_noad = ori_pref_noad(well_fit_bool);
ori_pref_noad_boot_avg = ori_pref_noad_boot_avg(well_fit_bool);
ori_pref_noad_boot_med = ori_pref_noad_boot_med(well_fit_bool);

plot(ori_pref_noad_boot_avg)
hold on
plot(ori_pref_noad_boot_med)
plot(ori_pref_noad)
legend('boot avg', 'boot med', 'fit')
