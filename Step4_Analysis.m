% Step 4: Analysis
clear;clc;
%% Setup Project Names
ProjectName = 'ThalamusSP_Project';
SUB = 'jcr70';
EXPER_SESSION = '20191016';
%% Add Pathways
BASEPATH = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\';
PROJECT  = [BASEPATH ProjectName '\'];

CODE_GENERAL = [BASEPATH 'Code\General\'];
%CODE_BRITTON = [BASEPATH 'Code\Packages\britton_code\code\'];
%CODE_BRITTON_PLOT = [BASEPATH 'Code\Packages\britton_code\other_code'];
DATA_STRUCTS = [PROJECT 'Data\jcr70\ephys\jcr70_20191016_920um_g0\'];

addpath(genpath(CODE_GENERAL));
addpath(genpath(DATA_STRUCTS));
%addpath(genpath(CODE_BRITTON));
%addpath(genpath(CODE_BRITTON_PLOT));

%% Load Structs
% Load data struct
load([DATA_STRUCTS 'data_struct_jcr70_20191016_920um_g0_t0.mat']);
load([DATA_STRUCTS 'lfp.mat']);
%% Setup
SetGraphDefaults;

% Set Color Palette
ColorPlt = [0,66,157;...
           46,89,168;...
           71,113,178;...
           93,138,189;...
           115,162,198;...
           138,188,207;...
           165,213,216;...
           197,237,223;...
           255,255,224];
ColorPlt = ColorPlt/255;

%% Rename variables
dig2.x = data_struct.trajectories.dig2_x;
dig2.y = data_struct.trajectories.dig2_y;
dig2.z = data_struct.trajectories.dig2_z;

dig4.x = data_struct.trajectories.dig4_x;
dig4.y = data_struct.trajectories.dig4_y;
dig4.z = data_struct.trajectories.dig4_z;

wrist.x = data_struct.trajectories.wrist_x;
wrist.y = data_struct.trajectories.wrist_y;
wrist.z = data_struct.trajectories.wrist_z;

pellet.x = data_struct.trajectories.pellet_x;
pellet.y = data_struct.trajectories.pellet_y;
pellet.z = data_struct.trajectories.pellet_z;

liftFrames = data_struct.jaaba.lift;
handOpenFrames = data_struct.jaaba.handopen;
grabFrames = data_struct.jaaba.grasp;
atMouthFrames = data_struct.jaaba.atmouth;
%% Trajectories
% Correct for pellet position
% add what frame to center around

%%%  Plot 1D trajectories for X, Y, Z  %%%
% Lift frames are plotted before 400 frames and after 1000 frames
figure;
sgtitle("1D Trajectories: Lift Aligned");
% Plot X
subplot(1,3,1);
plot1DTrajectories_fromXYZ(dig2.x, liftFrames);
title("X Position")
% Plot Y
subplot(1,3,2);
plot1DTrajectories_fromXYZ(dig2.y, liftFrames) %,'SampleSize', 50);
title("Y Position");
% Plot Z
subplot(1,3,3);
plot1DTrajectories_fromXYZ(dig2.z, liftFrames);%,'SampleSize',50);
title("Z Position");
%%%   Plot 3D trajectories   %%% WORK HERE _ HOW TO GET POINT AVG REACH
% LiftFrames will be the first plot point
figure;
subplot(1,2,1);
sgtitle("Reach Trajectory (Digit 4)")
% Plot Digit 4
plot3DTrajectories_fromXYZ(dig4.x,dig4.y,dig4.z,pellet,liftFrames,handOpenFrames,grabFrames,atMouthFrames,'SampleSize',50);
subplot(1,2,2);
% Plot Same thing without samples and without box or lines
plot3DTrajectories_fromXYZ(dig4.x,dig4.y,dig4.z,pellet,liftFrames,handOpenFrames,grabFrames,atMouthFrames,'SampleSize',50,'PlotMeanOnly',1);
% Another plot but make it a video that rotate
figure;
plot3DTrajectories_fromXYZ(dig4.x,dig4.y,dig4.z,pellet,liftFrames,handOpenFrames,grabFrames,atMouthFrames,'SampleSize',50,'PlotMeanOnly',1);
daspect([1,1,.3]);axis tight;
OptionZ.FrameRate=15;OptionZ.Duration=30;OptionZ.Periodic=true;
FileNameVideo = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\ThalamusSP_Project\Figures\jcr70\Behavior\3D_TrajectoryVideo.mp4';
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],FileNameVideo,OptionZ);


