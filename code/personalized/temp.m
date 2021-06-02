t1 = load('corr_well_fit_HVA.mat');
t2 = load('corr_well_fit_HVA_rerun.mat');

t1 = t1.ori_perc_all;
t2 = t2.ori_perc_all;

sum(t1==t2)