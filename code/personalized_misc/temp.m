t = squeeze(data_temp);
t_avg = mean(t(:,:,4000:4050), 3);
size(t_avg)

subplot(211)
imagesc(t(:,:,4050))
subplot(212)
imagesc(t_avg)