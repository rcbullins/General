function [scoreLabel] = specifyScoreLabel(score)
% PURPOSE
%   For specifing what trials to grab (success vs no success) and giving
%   each a name to be able to understand later. 

scoreLabel = num2str(score);

if strcmp(scoreLabel, '1')
    SCORE = 'idealSuccess';
elseif strcmp(scoreLabel, '0')
    SCORE = 'noSuccess';
elseif strcmp(scoreLabel, '-1')
    SCORE = 'noReach';
elseif strcmp(scoreLabel, '2')
    SCORE = 'eventualSuccess';
elseif strcmp(scoreLabel, 'all')
    score = [1 0 2];
    SCORE = 'allScores';
end

end