mask_np_overlay = sum(mask_np,3);
subplot(2,1,1)
imshow(mask_np_overlay)
subplot(2,1,2)
imshow(mask_cell)