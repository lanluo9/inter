function npSub_tc = get_cellpose_timecourse(data_reg, tif_name, LL_base, arg_date, imouse, run_str)

% t = load('cellpose_mask.mat');
t = load(tif_name);

mask_cell = t.cellpose_mask;

if max(mask_cell(:)) > 0
    
    while isempty(dir('*TCs_cellpose*.mat')) 
    % while cellpose timecourse does not exist, try generate it / wait for enough memory

        try
            
            mask_np = imCellNeuropil(mask_cell, 3, 5);
            
            % neuropil subtraction
            down = 5; % down sampling
            sz = size(data_reg);
            data_tc = stackGetTimeCourses(data_reg, mask_cell);
            disp('start stack_group_proj')
            data_reg_down  = stackGroupProject(data_reg,down);
            disp('finish stack_group_proj')
            data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
            nCells = size(data_tc,2)
            
            np_tc = zeros(sz(3),nCells);
            np_tc_down = zeros(floor(sz(3)./down), nCells);
            for i = 1:nCells
                np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
                np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
 			    log_every_n = floor(nCells ./ 10); % log less frequently
			    if mod(i, log_every_n)==0 
				    fprintf(['Cell #' num2str(i) '%s /n'])
                    disp(' ')
			    end
            end
            
            %get weights by maximizing skew
            ii= 0.01:0.01:1;
            x = zeros(length(ii), nCells);
            for i = 1:100
                x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
            end
            [max_skew, ind] =  max(x,[],1);
            np_w = 0.01*ind;
            npSub_tc = data_tc - bsxfun(@times,tcRemoveDC(np_tc),np_w);
            clear data_reg data_reg_down
            
            date = arg_date;
            save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], ...
                [date '_' imouse '_' run_str], ...
                [date '_' imouse '_' run_str '_TCs_cellpose.mat']), ...
                'data_tc', 'np_tc', 'npSub_tc', 'mask_cell', 'mask_np')
            clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
            disp(['TC extraction complete for ', run_str])
        
        catch
            disp('out of memory atm, retrying in a bit')
            pause(60)
        end

    end

else 
    disp('skipped cellpose TC extraction due to cellpose bug: ncell=0')
end

end