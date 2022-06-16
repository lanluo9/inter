%% search field name in struct

tt = fieldnames(input);
index = cellfun(@(x) any(contains(lower(x), 'grat')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input.', tt{id(i)}])
    fprintf('\n')
end

% % spatial freq
% input.tGratingSpatialFreqCPD
% input.gratingSpatialFreqCPD
% input.block2GratingSpatialFreqCPD
% input.photoMaskSpatialFreq

% % phase
% input.gratingStartingPhaseDeg