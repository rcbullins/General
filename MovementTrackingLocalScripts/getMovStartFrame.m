function [movStartFrames] = getMovStartFrame(trialScoreThr,traj,trialScores,trialBlock, varargin)
% PURPOSE
%   get start frame indexes of move time intiation.
% INPUTS
%   trialScoreThr (which trials you want to look at)
%       1 = success and one reach (ideal scenario)
%       2 = multiple reaches and success
%       0 = reach and grab without success
%       -1 = nothing
%       g = grooming
%   trialScores
%       vecotr of the trial score for each trial
%   trialBlock
%       trialIdxs.nbase (baseline trials only)
%       trialIdxs.npert (stimulation trials only)
%       trialIdxs.nwash (post stim trials only)
%   traj
%       trajectory output from get_traj_3D_v2
%   plotSanity (0 or 1) default = 1
% OUTPUS
%   movStartFrames
%       vector of all start frames for all trials
%   sanity check plot
%       plots sample trajectories with identified start time
% EXAMPLE
%   getMovStartFrame(1,traj, trialIdxs.trialScore,trialIdxs.nbase) will find
%   the start frames for all trials with a score of 1 during baseline
% HISTORY
%   Reagan Bullins 11.12.2021
%% Input Parser
p = inputParser;
addParameter(p,'plotSanity',0,@isnumeric);
addParameter(p,'allTrials',0,@isnumeric); % just find all trials, forget rest
parse(p,varargin{:});
plotSanity = p.Results.plotSanity;
find_allTrials = p.Results.allTrials;
%%
if find_allTrials == 0
    % Find the trials that have the specified trial score
    idx_specified = [nan];
    for itrialScore = 1:length(trialScoreThr)
        thisScore_idx = find(trialScores == trialScoreThr(itrialScore));
        idx_specified = [idx_specified thisScore_idx];
    end
    if length(trialScoreThr) > 1
        idx_specified(1) = [];
    end
    % Find where the block condition is met (baseline vs stimulation)
    [~,block_score_idx,~] = intersect(trialBlock,idx_specified);
    movStartFrames = NaN(1,length(block_score_idx));
    if plotSanity
        figure;
    end
    for itrial = 1:length(block_score_idx)
        thisTraj = traj(block_score_idx(itrial),:,:);
        z = squeeze(traj(block_score_idx(itrial),3,:));
        z_smooth = smoothdata(z,'movmedian',60);
        z_deriv = diff(z_smooth);
        max_z_deriv = max(z_deriv);
        thr = max_z_deriv*.5;
        overthr_idx = find(z_deriv >= thr);
        %             overbaseThr_idx = find(z < -30);
        %             [~,thr_idx,~] = intersect(overthr_idx,overbaseThr_idx);
        movStartFrames(1,itrial) = overthr_idx(1);
        if plotSanity && itrial <=16
            subplot(4,4,itrial)
            plot(z, 'Color',[.5 0 .5]);
            hold on;
            xline(movStartFrames(1,itrial));
        end
    end
else
     for itrial = 1:size(traj,1)
        thisTraj = traj(itrial,:,:);
        z = squeeze(traj(itrial,3,:));
        z_smooth = smoothdata(z,'movmedian',60);
        z_deriv = diff(z_smooth);
        max_z_deriv = max(z_deriv);
        thr = max_z_deriv*.5;
        overthr_idx = find(z_deriv >= thr);
        %             overbaseThr_idx = find(z < -30);
        %             [~,thr_idx,~] = intersect(overthr_idx,overbaseThr_idx);
        movStartFrames(1,itrial) = overthr_idx(1);
        if plotSanity && itrial <=16
            subplot(4,4,itrial)
            plot(z, 'Color',[.5 0 .5]);
            hold on;
            xline(movStartFrames(1,itrial));
        end
    end
end