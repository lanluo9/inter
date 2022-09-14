function [] = save_mat_as_tif(data_dfof_max)

stim_resp_gauss = data_dfof_max'; % gauss smooth each stim resp, take max
disp('saving dfof-max-gauss tif')

tif_file = 'cellpose_stim_resp_gauss.tif';
fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = 1;
B = numel(stim_resp_gauss)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(stim_resp_gauss(:,:,ct));
end
fTIF.close;
fprintf(1,repmat('\b',[1,msg]));
t=toc;
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

end