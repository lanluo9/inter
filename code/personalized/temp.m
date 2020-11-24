
x = []; y = [];
for ibin = 1 : nbin
    for iisi = 1 : nisi
        x(ibin) = bin_list(ibin);
        y(ibin, iisi) = mean(dis_pref_change(dis_pref_bin == x, iisi));
    end
end