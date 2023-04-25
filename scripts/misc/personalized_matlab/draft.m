icell = 1
icond = 3 % isi = 250 ms

data % resp to 8 ori 
ori_pref_cond(icell) % pref ori via curve fit


plot([0:22.5:180], [data, data(1)])
hold on
xline(ori_pref_cond(icell))