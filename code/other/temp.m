cd 200728_i1324_V1
mat_files = dir('*.mat'); 
for q = 1:length(mat_files) 
    load(mat_files(q).name); 
end