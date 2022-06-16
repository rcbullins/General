function [] = VideoCombiningPrep(VIDEO_DIR)
% PURPOSE
%   Preps video organization to be ran on cluster for JAABA training later.
%   Will prep all sessions specified in Directory_Animals Script.
% INPUT
%   Video folder of all raw side and front videos for all trials
% TODO:
%   Test line 40-45: makes copy of video for you if there is not a
%   copy folder already made. (Previously I hand copied the raw video folders
%   as an extra back up - but the code should now do it for you if a copy does
%   not exist)
% HISTORY
%   1/12/2022: Reagan Bullins
%% Find the movement trajectory

% Identify path with trk files (estimated position from
% tracker)

VIDEO_DIR_COPY = [VIDEO_DIR(1:end-1) ' - Copy/'];

if ~exist(VIDEO_DIR_COPY) %have not tested this
    % if video copy does not exist, make a copy
    mkdir(VIDEO_DIR_COPY);
    % Make the new movie and trk file name and directories
    copyfile(VIDEO_DIR, VIDEO_DIR_COPY);
end
setUpDir_Part1(VIDEO_DIR_COPY,'frontside',true);
end
%StartJAABA;