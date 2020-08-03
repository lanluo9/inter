%% get all cell trial trace by dir & isi

target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal
targ_start_list = unique(targ_start);
ngap = length(targ_start_list);

trace_cond = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        id_targ = find(targ_start == targ_start_list(igap));
        
        idx = intersect(intersect(id_targ, id_delta), id_adapter); % with-ad, specific isi & ori
        range_trace = [1 : max(trial_len)]; 
        trace_cond{icell, idelta, igap} = squeeze(tc_trials(icell, idx, range_trace)); % [ntrial, trial_len]
    end
end
end

%%
trace_no_ad = cell(ncell, ndelta); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        id_targ = find(targ_start == targ_start_list(igap));
        idx = intersect(intersect(id_targ, id_delta), id_noadapter); % no-ad, merge isi & specific ori
        range_trace = [1 : max(trial_len)]; 
        trace_cond{icell, idelta, igap} = squeeze(tc_trials(icell, idx, range_trace)); 
    end
end
end

cd C:\Users\lan\Documents\repos\inter\code
save trace_by_cond.mat trace_cond trace_no_ad

