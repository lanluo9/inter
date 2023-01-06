function [arr_thresholded, id_retained] = threshold_percentile(arr_1D, discard_perc_low, discard_perc_high)

    % threshold a 1D array by its percentile
    % input: 
    %   1D array
    %   top percentile to discard
    %   bottom percentile to discard
    % output:
    %   1D array, thresholded to discard top and bottom percentiles
    %   index in 1D array retained after thresholding [bool]
    
    [nrow, ncol] = size(arr_1D);
    assert(nrow==1 || ncol==1);
    
    thres_low = prctile(arr_1D, discard_perc_low);
    thres_high = prctile(arr_1D, 100 - discard_perc_high);
    
    id_retained = (arr_1D >= thres_low) & (arr_1D <= thres_high);
    n_discard = sum(~id_retained);
    n_total = length(id_retained);
    discarded_percent = round(n_discard / n_total * 100);
    fprintf('discarded %d percent of 1D array', discarded_percent)

    arr_thresholded = arr_1D(id_retained);

end