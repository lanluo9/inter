tt = fieldnames(input);
index = cellfun(@(x) any(contains(x, 'Ms')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input.', tt{id(i)}])
    fprintf('\n')
end


input.tStimOneGratingOnTimeMs
input.tStimTwoGratingOnTimeMs

input.mwStimOneOnMs
input.mwStimTwoOnMs
input.mwStimOneOffMs
input.mwStimTwoOffMs

input.stimOneGratingOnTimeMs
input.stimTwoGratingOnTimeMs
