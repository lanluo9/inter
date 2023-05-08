function ori_pref = fit_tuning_half_trials(dfof_align_tg, save_flag)

% fit tuning twice, on 2 halves of no adapter trials
% input: dfof_align_tg. output: ori_pref, where nrun = 2

global ncell nori id_ori ori_list id_isi3 range_base range_resp

nisi = length(id_isi3);
nrun = 2;

dfof_avg = pi * ones(ncell, nori, nisi); 
fit_param = pi * ones(ncell, 7, nisi); 
ori_pref = pi * ones(ncell, nisi);
ori_pref_runs = pi * ones(ncell, nisi, nrun);

for irun = 1 : nrun
for icell = 1 : ncell
for iisi =  1 : length(id_isi3) 
    for iori = 1 : nori
        idx = intersect(id_ori{iori}, id_isi3{iisi});
        ntrials_cond = length(idx);
        ntrial_draw = round(ntrials_cond * 0.5);
        if irun == 1
            idx_run = idx(1 : ntrial_draw);
        elseif irun == 2
            idx_run = idx(ntrial_draw : end);
        end
%         idx_run = randsample(idx, bootstrap_draw, 0); % without replacement
%         idx_run_complement = setdiff(idx, idx_run); % the other half of ntrials, no overlap

        base_win = squeeze(dfof_align_tg(icell, idx_run, range_base)); base_win = mean(base_win, 2);
        resp_win = squeeze(dfof_align_tg(icell, idx_run, range_resp)); resp_win = mean(resp_win, 2);
        dfof_avg(icell, iori, iisi) = mean( resp_win - base_win );
    end

    data = dfof_avg(icell, :, iisi); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(deg2rad(ori_list), data);
    fit_param(icell, :, iisi) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%   icell, baseline offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2

    ori_pref = rad2deg(u1_hat);
    ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
    ori_pref_runs(icell, iisi, irun) = ori_pref;
end
end
end

if save_flag
    save fit_tuning_half_trials.mat ori_pref_runs; 
end


