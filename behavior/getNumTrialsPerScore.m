function [scoreCt] = getNumTrialsPerScore(trialIdxs, experCondition)
% PURPOSE
%   Get how many trials were for each score (success vs no success).
% INPUT
%   trialIdxs: the index of each trial with its score on accuracy.
%   experCondition: stim/washout/base/etc. indices of these trials.
% OUTPUT
%   Count of trials per score.
scoreCt.idealSuccess = length(find(trialIdxs.trialScore(experCondition.trialIdxs) == 1));
scoreCt.eventualSuccess = length(find(trialIdxs.trialScore(experCondition.trialIdxs) == 2));
scoreCt.noSuccess = length(find(trialIdxs.trialScore(experCondition.trialIdxs) == 0));
scoreCt.noReach = length(find(trialIdxs.trialScore(experCondition.trialIdxs) == -1));
scoreCt.grooming = length(find(trialIdxs.trialScore(experCondition.trialIdxs) == 'g'));
end