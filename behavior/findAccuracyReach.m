function [acc_performance] = findAccuracyReach(numPerformance, varargin)
% PURPOSE
%   Find Reach performance by calculating accuracy of reach task. Create
%   plot if plotAccuracy = 1.
% INPUT
%   numPerformance (Struct)         num trials per category.
%           .idealSuccess (Int)
%           .eventualSuccess
%           .noReach
% OUTPUT
%    Accuracy per category and bar plot.
%% Input Parsers
p = inputParser;
addParameter(p,'graphTitle','Performance Accuracy',@isstr);
addParameter(p,'plotAccuracy',1,@isnumeric);

parse(p,varargin{:});
graphTitle         = p.Results.graphTitle;
plotAccuracy        = p.Results.plotAccuracy;
%% 
numTotal = numPerformance.idealSuccess+numPerformance.eventualSuccess+numPerformance.noSuccess+numPerformance.noReach;
numReaching = numPerformance.idealSuccess+numPerformance.eventualSuccess+numPerformance.noSuccess;
acc_performance = [numPerformance.idealSuccess/numReaching; numPerformance.eventualSuccess/numReaching ;numPerformance.noSuccess/numReaching ;numPerformance.noReach/numTotal];

if plotAccuracy == 1
    figure;
    h = bar(acc_performance*100,'FaceColor','flat');
    colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
    h.CData(1,:) = colorScheme(1,:);
    h.CData(2,:) = colorScheme(2,:);
    h.CData(3,:) = colorScheme(3,:);
    h.CData(4,:) = colorScheme(4,:);
    l = cell(1,4);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach';
    set(gca,'xticklabel', l);
    ylabel('Accuracy (%)');
    title(graphTitle);
end