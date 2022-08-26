function [] = plot1DTrajectories_fromXYZ(pos, liftFrames,varargin)
%% Behavioral - Motor reaching kinematics analysis script
% PURPOSE
%   Plot 1D trajectories of animal movement
% INPUT
%   PlotFrameBefore  how many frames to plot before lift (default 400)
%   PlotFrameAfter   how many frames to plot after lift (default 1000)
%   pos        position information, x,y, or z information over time
%   liftFrames  which frame lift occured in for each trial (1 x num trial)
% OUTPUT
%   sample trial trajectory information 
% HISTORY
%   Reagan Bullins: 8.24.2022
%% Input parsers
p = inputParser;
addParameter(p,'Color', [1, 0, 0], @isvector)
addParameter(p,'SampleSize',size(pos,1), @isnumeric);
addParameter(p,'PlotFrameBefore',400, @isnumeric);
addParameter(p,'PlotFrameAfter',1000, @isnumeric);
parse(p,varargin{:});
Color = p.Results.Color;
SampleSize = p.Results.SampleSize;
PlotFrameBefore = p.Results.PlotFrameBefore;
PlotFrameAfter = p.Results.PlotFrameAfter;
%% Plot individual trials on the specified axis
%plot trials (sample of trials)
tmp_mean = [];

for i=1:size(pos,1)
    %if no lift for this trial, then skip it
    if isnan(liftFrames(1,i))
        continue;
    else
    % else plot the lift, however many frames before and after lift
        if i <= SampleSize
            plot((-PlotFrameBefore:PlotFrameAfter)*2,pos(i, liftFrames(1,i)-PlotFrameBefore:liftFrames(1,i)+PlotFrameAfter), 'Color', [105/255 105/255 105/255 .1]);
        end
        tmp_mean = [tmp_mean; pos(i, liftFrames(1,i)-PlotFrameBefore:liftFrames(1,i)+PlotFrameAfter)];
    end
    hold on;
end
xlabel('Time (ms)');

plot((-PlotFrameBefore:PlotFrameAfter)*2,mean(tmp_mean),'r')

% Find where lift does not happen and remove from mean
%non_nan_idx = find(isnan(liftFrames) == 0);
% plot mean of all samples, not just sample size
%pos_mean = mean(pos(non_nan_idx',(liftFrames(1,non_nan_idx)-PlotFrameBefore):(liftFrames(1,non_nan_idx)+PlotFrameAfter)));
%plot((-PlotFrameBefore:PlotFrameAfter)*2,pos_mean,'Color', Color);
%hold on;


axis square;
end