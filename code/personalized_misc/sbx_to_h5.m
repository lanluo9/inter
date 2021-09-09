function sbx_to_h5(fname,varargin)
 
% sbx2h5: Generates h5 file from sbx
% input: fname is a string, i.e. '002_000_000' for 002_000_000.sbx
% https://scanbox.org/2018/08/29/using-suite2p-with-scanbox/

% cd Z:\All_Staff\home\lan\Data\2P_images\i1324\200728\003
% sbx_to_h5('003_000_000')

fnh = [fname ,'.h5']; 
z = sbxread(fname,1,1);
global info;
 
if(nargin>1)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end
 
k = 0;
done = 0; 
blksize = 1000; % block size original 200 
to_read = min(blksize,N-k);
 
while(~done && to_read>0)
try
    q = sbxread(fname,k,to_read);
    q = squeeze(q(1,:,:,:)); % extract green channel only
    q = permute(q,[2 1 3]);
%     nline = 512;
    nline = 264;
    if(k==0) % assumes frames with 512 lines. w/ frame_rate = 30, should be 264 lines
%       30.04 * 264 ~ 15.49 * 512
%       data = 264 x 796 x 100000 uint16
        h5create(fnh,'/data',[796 nline Inf],'DataType','uint16','ChunkSize',[796 nline to_read]);
        h5write(fnh,'/data',q,[1 1 1],[796 nline to_read]);
        f = waitbar(0,'Converting to hdf5');
    else
        h5write(fnh,'/data',q,[1 1 k+1],[796 nline to_read]);
    end
catch
    done = 1;
    delete(f);
end

k = k + to_read;
to_read = min(blksize,N-k);
waitbar(k/N,f);
end
 
delete(f);
end