% Step 4: Analysis
% RB
clear;clc;
%% Setup Project Names
ProjectName = 'ThalamusSP_Project';
SUB = 'jcr70';
EXPER_SESSION = '20191016';
%% Add Pathways
BASEPATH = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\';
PROJECT  = [BASEPATH ProjectName '\'];

CODE_GENERAL = [BASEPATH 'Code\General\'];
DATA_STRUCTS = [PROJECT 'Data\jcr70\ephys\jcr70_20191016_920um_g0\'];

addpath(genpath(CODE_GENERAL));
addpath(genpath(DATA_STRUCTS));

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
% Behavior
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
% Neural
% load in all z score data for all trials
M1_z_score = {};
    for itrial = 1:length(data_struct.trials)
        M1_z_score{itrial} = data_struct.trials(itrial).M1_z_score;
    end
trial_time = data_struct.trials.time;
%% Reaching Trajectories
% Correct for pellet position - if multiple days...
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
%% Neural Trajectories
% Select trials with full labeled sequence
    non_nan_idx_lift = find(~isnan(liftFrames));
    non_nan_idx_hand = find(~isnan(handOpenFrames));
    non_nan_idx_grasp = find(~isnan(grabFrames));
    non_nan_idx_mouth = find(~isnan(atMouthFrames));
    % Find all non nan idx for all behaviors
    non_nan_idx  = intersect(intersect(intersect(non_nan_idx_lift,non_nan_idx_hand),non_nan_idx_grasp),non_nan_idx_mouth);
    % Get the lift average frame
    lift_mean_frame = round(mean(liftFrames(1,non_nan_idx)));
    mouth_mean_frame = round(mean(atMouthFrames(1,non_nan_idx)));
    mouth_inRespectToLift = mouth_mean_frame-lift_mean_frame;
% Average over trials
    meanMatrix = zeros(size(M1_z_score{1},1),size(M1_z_score{1},2));
    for i = 1:size(M1_z_score,2)
        currentMatrix = M1_z_score{i}; 
        meanMatrix = (meanMatrix + currentMatrix);
    end
    avg_M1_z_score = meanMatrix/size(M1_z_score,2);
% PCA over z scored data, specific trials
[coeff,score,latent] = pca(avg_M1_z_score','NumComponents',3);
trialFrameStart = find(trial_time == 0);
reachFrames = lift_mean_frame+trialFrameStart:mouth_mean_frame+trialFrameStart;
behaviorFrames = trialFrameStart:length(trial_time);
PlotFrameAfter = 1000;
 
figure;
sgtitle('Neural Trajectory: PC Space')
subplot(1,2,1);
plot3(score(behaviorFrames,1),score(behaviorFrames,2),score(behaviorFrames,3));
hold on;
plot3(score(reachFrames,1),score(reachFrames,2),score(reachFrames,3),'r');
scatter3(score(reachFrames(1),1),score(reachFrames(1),2),score(reachFrames(1),3),'o','filled','MarkerEdgeColor',ColorPlt(1,:),'MarkerFaceColor', ColorPlt(1,:));
txt_lift = 'Lift';
text(score(reachFrames(1),1),score(reachFrames(1),2),score(reachFrames(1),3),txt_lift);
subplot(1,2,2);
plot3(score(reachFrames,1),score(reachFrames,2),score(reachFrames,3),'r');
hold on;
scatter3(score(reachFrames(1),1),score(reachFrames(1),2),score(reachFrames(1),3),'o','filled','MarkerEdgeColor',ColorPlt(1,:),'MarkerFaceColor', ColorPlt(1,:));
txt_lift = 'Lift';
text(score(reachFrames(1),1),score(reachFrames(1),2),score(reachFrames(1),3),txt_lift);
axis off;
box off;
%% Spiking Characteristics