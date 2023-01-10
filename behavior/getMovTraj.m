function [traj, traj_smooth] = getMovTraj(data_dir, CODE_CALIB,varargin)
% PURPOSE
%   If the movement trajectory file does not exist for this bodypart, make
%   the file. Load the trajectory and save.

% INPUTS
%   trajfile: path to trajectory mat file
%   TRK: path to trk file output from APT
%   CODE_CALIB: path to calibration file for cameras
%   NumVidFrames: how many frames to analyze in videos
%   digit2trk: body part to track (1=1st digit, 2=3rd digit,3=wrist,4=pellet)
        % NOTE: this could be different - depending on how you label in APT

% OUTPUT
%   traj: 3D matrix (trials x coordinate xyz x frames); 
%% Input Parsers
p = inputParser;
addParameter(p,'NumVidFrames',2500,@isnumeric);
addParameter(p,'dig2trk',1,@isnumeric);
addParameter(p,'smoothTrace',1,@isnumeric);
parse(p,varargin{:});

NumVidFrames              = p.Results.NumVidFrames;
dig2trk                   = p.Results.dig2trk;
%% If the file exists, load it, if it does not exist, create it
% if exist(trajFile,'file')==2
%     load(trajFile);
% else
    [traj, ~] = get_traj_3D_RB(data_dir,CODE_CALIB,NumVidFrames,dig2trk);
end
