load carbig
p = anovan(MPG,{org when},'model',2,'varnames',{'origin','mfg date'})


p_anova_n = anovan(adp_ratio_cell,{area mouse isi}, 'model',3, 'varnames',{'area', 'mouse', 'isi'})
