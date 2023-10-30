function z_make_fig6f_cleaned()
% This function generates Figure 6f based on data from multiple files.

%% Data Preparation
% List of file names
filename = {'170323_i689_runs-002-003', '170323_i696_runs-002-003', '170324_i674_runs-002-003', ...
            '170327_i684_runs-002-003', '170503_i711_runs-002-003', '170503_i712_runs-002-003', ...
            '170510_i574_runs-002-003', '170808_i720_runs-002-003', '170810_i738_runs-002-003', ...
            '170811_i739_runs-002-003', '170816_i745_runs-002-003', '170826_i746_runs-002-003'};
        
PCmax = 15;
train = 3;

% Initialize data structures
DVAll.dataset = [];
DVAll.Y = [];
DVAll.cond = [];
DVAll.PV = [];

% Loop through files
for n = 1:length(filename)
    % Load data from files
    load([filename{n}, '_newFits.mat']);
    
    [~, loc] = max(ori_fit);
    prefs{n} = loc / 180 * 2 * pi;
    idxn{n} = find(theta_90 < 22.5);
    
    % Display the number of good units using theta_90
    disp(['Dataset ', num2str(n), ' has ', num2str(length(idxn{n})), ' good units using theta_90!']);
    
    prefs{n} = prefs{n}(idxn{n});
    f{n} = ori_fit(:, idxn{n});
    kappa{n} = abs(fft(log(f{n})));
    kappa{n} = kappa{n}(end, :);

    % Data preparation for analysis
    for j = 1:max(train, 2)
        data{j}.X = [];
        data{j}.Y = [];

        for k = 1:8
            data{j}.X = [data{j}.X; ppResp{j, k}'];
            data{j}.Y = [data{j}.Y; k * ones(size(ppResp{j, k}, 2), 1)];
        end

        % Remove rows with NaN values
        idx = ~any(isnan(data{j}.X), 2);
        data{j}.X = data{j}.X(idx, :);
        data{j}.Y = data{j}.Y(idx, 1);

        data{j}.X = data{j}.X(:, idxn{n});
        data{j}.Xraw = [data{j}.X];
    end

    % Calculate principal components
    CC = cov([data{1}.X; data{2}.X]);
    [V, D] = eig(CC, 'vector');
    PCs = min(PCmax, size(data{j}.X, 2));
    D = D(max(size(data{1}.X, 2) - PCs + 1, 1):end);
    V = V(:, max(size(data{1}.X, 2) - PCs + 1, 1):end);

    dataAll.X = [];
    dataAll.Xraw = [];
    dataAll.Y = [];
    dataAll.cond = [];
    for j = 1:max(train, 2)
        data{j}.X = data{j}.X * V * diag(1 ./ sqrt(D));
        data{j}.X = data{j}.X(:, max(size(data{j}.X, 2) - PCs + 1, 1):end);
        data{j}.X = [data{j}.X, ones(size(data{j}.X, 1), 1)];
        dataAll.X = [dataAll.X; data{j}.X];
        dataAll.Xraw = [dataAll.Xraw; data{j}.Xraw];
        dataAll.Y = [dataAll.Y; data{j}.Y];
        dataAll.cond = [dataAll.cond; j * ones(size(data{j}.Y))];
    end

    % Calculate DVs
    [~, DVAlltemp, Btemp, ~] = getDVs(data, dataAll, 0, prefs{n}(:,:), kappa{n}(:,:), f{n}(:,:));

    for j = 1:max(train, 2)
        Btemp{j} = Btemp{j}(1:end - 1, 1);
    end

    % Remove bias from estimators
    for j = 1:max(train, 2)
        DVAll.dataset = [DVAll.dataset; n * ones(size(DVAlltemp{j}.opt))];
        DVAll.Y = [DVAll.Y; data{j}.Y];
        DVAll.cond = [DVAll.cond; j * ones(size(data{j}.Y));];
        DVAll.PV = [DVAll.PV; DVAlltemp{j}.PV];
    end
end

%% Plotting
not8 = [1, 2, 3, 5, 6, 7];
NDC = 500;
dv = [0:NDC] / NDC;
usedatasets = [2:4, 6:12];
kk = 0;

for dataset = usedatasets
    kk = kk + 1;
    for j = 1:2
        for difficulty = 1:5
            idxfp = (DVAll.Y == 8 & DVAll.cond == 1 & logical(sum(DVAll.dataset == dataset, 2)));
            switch difficulty
                case 1
                    not8 = 8;
                case 2
                    not8 = [1, 7];
                case 3
                    not8 = [2, 6];
                case 4
                    not8 = [3, 5];
                case 5
                    not8 = 4;
            end
            idxcd = (logical(sum(DVAll.Y == not8, 2)) & DVAll.cond == j & logical(sum(DVAll.dataset == dataset, 2)));
            DVtemp = abs(getfield(DVAll, 'PV'));
            CD = mean(DVtemp(idxcd) > dv);
            FP = mean(DVtemp(idxfp) >= dv);

            AUROC{k}(kk, j, difficulty) = -trapz(FP, CD);
        end
    end
end

xa = [0, 22.5, 45, 67.5, 90];

figure
errorbar(xa, mean(squeeze(AUROC{k}(:, 1, :))), std(squeeze(AUROC{k}(:, 1, :))) / sqrt(length(usedatasets)), 'b')
hold on
errorbar(xa, mean(squeeze(AUROC{k}(:, 2, :))), std(squeeze(AUROC{k}(:, 2, :))) / sqrt(length(usedatasets)), 'r')
hold off
title('PV')
ylabel('AUROC')
xlabel('Orientation difference')
axis([0, 90, 0.4, 1])
legend('250', '750')
end
