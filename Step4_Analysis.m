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
% Load in all spikes 
spikes = {};
    for itrial = 1:length(data_struct.trials)
        spikes{itrial} = data_struct.trials(itrial).M1_tstamps;
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
box off;
axis off;
% Plot Y
subplot(1,3,2);
plot1DTrajectories_fromXYZ(dig2.y, liftFrames) %,'SampleSize', 50);
title("Y Position");
box off;
axis off;
% Plot Z
subplot(1,3,3);
plot1DTrajectories_fromXYZ(dig2.z, liftFrames);%,'SampleSize',50);
title("Z Position");
box off;
axis off;
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
% Plot PCA in 3 dimensions and label start of lift movement
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
% Plot first three PCA each by themselves over time (plot lift)
figure;
subplot(1,3,1);
plot(score(behaviorFrames,1));
hold on;
x_pt = find(score(behaviorFrames,1) == score(reachFrames(1),1));
scatter(x_pt, score(reachFrames(1),1),'o','filled','MarkerEdgeColor',ColorPlt(1,:),'MarkerFaceColor', ColorPlt(1,:));
title("PC 1");
box off;
subplot(1,3,2);
plot(score(behaviorFrames,2))
hold on;
x_pt = find(score(behaviorFrames,2) == score(reachFrames(1),2));
scatter(x_pt, score(reachFrames(1),2),'o','filled','MarkerEdgeColor',ColorPlt(1,:),'MarkerFaceColor', ColorPlt(1,:));
title("PC 2");
box off;
subplot(1,3,3);
plot(score(behaviorFrames,3));
hold on;
x_pt = find(score(behaviorFrames,3) == score(reachFrames(1),3));
scatter(x_pt, score(reachFrames(1),3),'o','filled','MarkerEdgeColor',ColorPlt(1,:),'MarkerFaceColor', ColorPlt(1,:));
title("PC 3");
box off;
%% Spiking Characteristics
numFramesPlot = 1501;
num_neurons = size(avg_M1_z_score,1);
% heat plot: each row a neuron x time (average z score)
figure;
% lift align all trials and all neurons
lift_align_neurons = {};
%for each trial
for itrial = 1:size(M1_z_score,2)
    this_trial_neurons = M1_z_score{itrial};
    if isnan(liftFrames(1,itrial))== 1
        lift_align_neurons{itrial} = [];
        continue;
    end
    this_aligned = this_trial_neurons(:, trialFrameStart+liftFrames(1,itrial) - 500:trialFrameStart+liftFrames(1,itrial)+1000);
    lift_align_neurons{itrial} = this_aligned;
end

% Average over all trials
lift_align_activity = zeros(num_neurons,numFramesPlot);
for ineuron = 1:num_neurons
    neuron_all = zeros(size(lift_align_neurons,2),numFramesPlot);
    for itrial = 1:size(M1_z_score,2)
        if isempty(lift_align_neurons{itrial})== 1
            neuron_all(itrial,:) = zeros(1,numFramesPlot);
            continue
        end
        neuron_all(itrial,:) = lift_align_neurons{itrial}(ineuron,:);
    end
    lift_align_activity(ineuron,:) = mean(neuron_all);
end

% Make graph
figure;
imagesc(lift_align_activity)

%% heatmap Raster part 2, not z scored

% lift align all trials and all neurons
lift_align_neurons = {};
%for each trial
for itrial = 1:size(M1_z_score,2)
    this_trial_neurons = spikes{itrial};
    if isnan(liftFrames(1,itrial))== 1
        lift_align_neurons{itrial} = [];
        continue;
    end
    this_aligned = this_trial_neurons(:, trialFrameStart+liftFrames(1,itrial) - 500:trialFrameStart+liftFrames(1,itrial)+1000);
    lift_align_neurons{itrial} = this_aligned;
end

% Average over all trials
lift_align_activity = zeros(num_neurons,numFramesPlot);
for ineuron = 1:num_neurons
    neuron_all = zeros(size(lift_align_neurons,2),numFramesPlot);
    for itrial = 1:size(M1_z_score,2)
        if isempty(lift_align_neurons{itrial})== 1
            neuron_all(itrial,:) = zeros(1,numFramesPlot);
            continue
        end
        neuron_all(itrial,:) = lift_align_neurons{itrial}(ineuron,:);
    end
    lift_align_activity(ineuron,:) = mean(neuron_all);
end

% Make graph
figure;
imagesc(lift_align_activity)
%% Raster with avg over top

figure;
for ineuron = 4%
    %subplot(5,8,ineuron)
    for itrial = 1:size(spikes,2)
        idx_spikes = find(spikes{itrial}(ineuron,:)>=1);

        scatter(trial_time(idx_spikes),ones(1,length(idx_spikes))*itrial,'.','Color','k')
        hold on
    end
end

