% visdiff('quickRet_LL.m','quickRet.m')
% visdiff('behavConstsAV.m','behavConstsAV_LL.m')

load('200720_1323_runs-003_reg_shifts.mat')

input.tBlock2TrialNumber % no-adapter trials
input.tGratingDirectionDeg  % dir
input.gratingDirectionDeg  %
input.tTotalStimFrames  % adapter vs no-adapter trial nframe
input.tFramesOff  %??

input.itiTimeMs % trial len
input.stimOnTimeMs
input.stimOffTimeMs
input.stimOffTimeIntervalMs
input.targetOnTimeMs

