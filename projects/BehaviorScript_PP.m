%% Behavior Thalamus Single Pulse
%% Clean workspace
clear;
close ALL;
clc;
%%%
%Behavior Template Script
% PURPOSE
%   Run all necessary behavior analysis for dataset. This
%   script will create 2D and 3D trajectories of reach to grab movements,
%   calculate accuracy, look at endpoints. Will analyze all sessions
%   specified in Directory_Animals Script.
% DEPENDENCIES
%   Britton Code
%       Animal Tracking Calibration Files
%       Original get_traj_3D code
%   Trk files
%       Output from APT tracking (estimated x,y,z position of body parts)
%           - One trk file per trial (reach to grab trial)
%   JAABA output
%       JAABA is a more accurate way to define parts of the reach to grab
%       movement. Start with videoCombiningPrep --> run over cluster -->
%       train on JAABA (Matlab 2014 StartJAABA in JAABAST toolbox)
%   Excel experimental session information sheets
%       One for each subject
%       Created at the time of experiment.
%       Each session: clear background = baseline, colored background = stim
%       General session sequence: baseline > stim VA ChrRhp > washout
%   Directory_Animals.m
%       Script outlining experimental info: sessions and subjects to run
%       analysis over. Consult excel sheets for identifying trials as
%       baseline, stimulation, or washout.
%   Reagan_code for harmaline project
%       required functions to run script
%       getStartFrame (of movement)
%       getEndFrame (of movement)
%       plot1DTrajectories
%       plot3DTrajectories
%       plotOverlay1D
%       plotOverlay3D
%       Behavioral_Comparison_Conditions (plots 3D trajectories comparing conditions
%       plots endpoints only, aligns
%       endpoints based on pellet location)
%       Behavioral_Metrics (plots accuracy comparing conditions)
% OUTPUTS
%       See section below to specify what graphs to output: Samples
%           plotSampleTrajectories (plot 2D and 3D sample trajectories)
%           plot1DOverlayTrajectories (plot all 2D trajectories on same
%               graph)
%           plot3DOverlayTrajectories (plot all 3D trajectories on same
%               graph)
%           plotAccuracy(plot individual accuracy for each condition)
%       Comparison graphs (will always make these)
%           Comparing accuracy on reach movement
%           Comparing 3D trajectories (from movement initation to pellet
%               grab)
%               - Another plotting endpoints only. Dots are pellet locations
%               - Another aligning all condition sessions using the
%                    location of the pellet. Then plotting endpoints as success
%                    (o) or failures (x)
% TO DO
%       Change JAABA output file load in (it does it the old way- not the
%       new way)
% HISTORY
%   11/23/2021 Reagan Bullins
%% Experiment
SUB = 'rb2';
EXPER_SESSION = '20221208';

%% Set Paths
BASEPATH = 'C:/Users/bullinsr/OneDrive - University of North Carolina at Chapel Hill/Hantman_Lab/';
CODE_REAGAN = [BASEPATH 'Code/General/'];
RAWDATA_BASEPATH = 'D:/rbullins/Data/';

SPIKE_GLX = [BASEPATH '/Code/Packages/SpikeGLX_Datafile_Tools/']
CODE_CALIB = ['D:/rbullins/Data/' SUB '/Calib_Results_stereo.mat'];
% CODE_BRITTON = [RAWDATA_BASEPATH 'Code/britton_code/code/matlab_britton/'];
% CODE_BRITTON_PLOT = [RAWDATA_BASEPATH 'Code/britton_code/other_code'];
% CODE_PROCESS_EVENTS = [RAWDATA_BASEPATH 'Code/process_events/'];

JAABA_OUTPUT = [RAWDATA_BASEPATH SUB '/' EXPER_SESSION '/JAABA_Combo/'];

addpath(genpath(JAABA_OUTPUT));
addpath(genpath(CODE_REAGAN));
addpath(genpath(SPIKE_GLX))
% addpath(genpath(CODE_BRITTON));
% addpath(genpath(CODE_BRITTON_PLOT));
% addpath(genpath(CODE_PROCESS_EVENTS));
%% Score - what trials to grab
score = 'all';
scoreLabel = specifyScoreLabel(score);

%% Find the movement trajectory

% Identify path with trk files (estimated position from
% tracker)
TRK = [RAWDATA_BASEPATH SUB '/' EXPER_SESSION '/trk/'];
addpath(TRK);
%% Load position data
% For digit 1
% trajFile1 = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories1.mat']);
traj = getMovTraj(TRK, CODE_CALIB,'dig2trk',1); %(trajFile1, TRK, CODE_CALIB,'digit2trk',1);

% For digit 2
% trajFile2 = ([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_trajectories2.mat']);
traj2 = getMovTraj(TRK, CODE_CALIB,'dig2trk',2);
%% Find any missing files
videos = [RAWDATA_BASEPATH SUB '/' EXPER_SESSION '/videos/'];
all_vids = dir([videos,'*front*.avi']);
vid_sec = zeros(1,length(all_vids));
sec_diff = zeros(1,length(all_vids)-1);
for ivideo = 1:length(all_vids)
   FileInfo = dir([videos all_vids(ivideo).name]);
   [year,month,day,hr,min,sec] = datevec(FileInfo.datenum);
   vid_sec(1,ivideo) = sec;
   if ivideo > 1
       if vid_sec(1,ivideo) < vid_sec(1,ivideo-1) %hit new minute
           sec_diff(1,ivideo-1) = vid_sec(1,ivideo) + (60 - vid_sec(1,ivideo-1));
       else
           sec_diff(1,ivideo-1) = vid_sec(1,ivideo) - vid_sec(1,ivideo-1);
       end
   end
end
suspect_indicies = find(sec_diff > 25); %plus one to this index for where the missing video should be
missing_vid = 66; % remove trial 66 in neural & dac stuff
%% Set trials
% Define Trials by baseline, stim, washout

% binary file
binName = [SUB '_' EXPER_SESSION '_g0_t0.nidq.bin'];
path = [RAWDATA_BASEPATH SUB '/' EXPER_SESSION '/' SUB '_' EXPER_SESSION '_g0/'];

% Parse the corresponding metafile
meta = ReadMeta(binName, path);

% Get all data
nSamp = str2num(meta.fileTimeSecs) * SampRate(meta);

dataArray = ReadBin(0, nSamp, meta, binName, path);

% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = 1;

if strcmp(meta.typeThis, 'imec')
    dataArray = GainCorrectIM(dataArray, [ch], meta);
else
    dataArray = GainCorrectNI(dataArray, [ch], meta);
end

% 7 Channels---- 3: Cue, 4: Camera, 6: Laser ----
cam_diff = diff(dataArray(4,:));

cam_diff(cam_diff<1000) = 0;
cam_diff(cam_diff>=1000)= 1;

smooth_data = movmean(cam_diff,30);
smooth_data(smooth_data > 0) = 1;

cam_ss = diff(smooth_data);

cam_start = find(cam_ss == 1) + 1;
cam_stop = find(cam_ss == -1);
cam_idx(1,:) = cam_start;
cam_idx(2,:) = cam_stop;

% Get all the cue times
cue = dataArray(3,:);
cue(cue <= 6000) = 0;
cue(cue > 6000) = 1;
cue_diff = diff(cue);

cue_start = find(cam_diff == 1) + 1;

% Chop up array via start and stop of camera
laser = dataArray(6,:);
numTrials = size(cam_idx,2);
if numTrials ~= size(cue_start,2)
    disp("Number of trials different for cue and camera")
end

for itrial = 1:numTrials
    laser_trial{itrial} = laser(1,cam_idx(1,itrial):cam_idx(2,itrial));
end
bt = 1;
p_20 = 1;
p_40 = 1;
s_20 = 1;
s_40 = 1;
s_10 = 1;
for itrial = 1:numTrials
    %%HERE WORKING
    if sum(b == a) == 1 %if this is a missing video, then skip it
        continue;
    end
    if max(laser_trial{itrial}) < 1000
        ctrl_trials{bt} = itrial;
        bt = bt+1;
    else
        if meanfreq(laser_trial{itrial}) >=.0039 && meanfreq(laser_trial{itrial}) <= .0043
            % 20 Hz stimulation Sine Expression
            sine_20{s_20} = itrial;
            s_20 = s_20+1;
        elseif meanfreq(laser_trial{itrial}) >=.0079 && meanfreq(laser_trial{itrial}) <= .0082
            % 40 Hz stimulation Sine Expression
            sine_40{s_40} = itrial;
            s_40 = s_40+1;
        elseif meanfreq(laser_trial{itrial}) >=.27 && meanfreq(laser_trial{itrial}) <= .3
            % 20 Hz stimulation Pulse
            pulse_20{p_20} = itrial;
            p_20 = p_20+1;
        elseif meanfreq(laser_trial{itrial}) >=.037 && meanfreq(laser_trial{itrial}) <= .04
            % 40 Hz stimulation Pulse
            pulse_40{p_40} = itrial;
            p_40 = p_40+1;
        elseif meanfreq(laser_trial{itrial}) >=.002 && meanfreq(laser_trial{itrial}) <= .0023
            sine_10{s_10} = itrial;
            s_10 = s_10 + 1;
        else
            disp(["Trial " num2str(itrial) ": Not Classified"]);
        end
    end
end
clear s_10 s_20 s_40 p_20 p_40 cam_start cam_stop cue cam_ss cam_diff cue_diff cue_times laser smooth_data cam_start_real binName
%% Set Graph Defaults - I like arial font :)
SetGraphDefaults;
%% Find accuracy for baseline and save mat
% Find total number of baseline trials
% numTotalBL = length(base.trialIdxs);
numTotalBL = numBL.idealSuccess+numBL.eventualSuccess+numBL.noSuccess+numBL.noReach;
numReaching = numBL.idealSuccess+numBL.eventualSuccess+numBL.noSuccess;
% Find number of each trial type
y_BL_acc = [numBL.idealSuccess/numReaching; numBL.eventualSuccess/numReaching;...
    numBL.noSuccess/numReaching ;numBL.noReach/numTotalBL];
% Save the accuracy as a mat file to be called upon later in
% Behavioral_Metrics
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baselineAccuracy.mat'],'y_BL_acc');
% Plot Accuracy (will be plotted in conjunction with other
% experimental conditions later)
if plotAccuracy == 1
    h = bar(y_BL_acc*100,'FaceColor','flat');
    colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
    h.CData(1,:) = colorScheme(1,:);
    h.CData(2,:) = colorScheme(2,:);
    h.CData(3,:) = colorScheme(3,:);
    h.CData(4,:) = colorScheme(4,:);
    l = cell(1,4);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach';
    set(gca,'xticklabel', l);
    ylabel('Accuracy (%)');
    title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' :BL Accuracy']);
end
%% Find stimulation accruacy and save mat
numTotalPert = numPert.idealSuccess+numPert.eventualSuccess+numPert.noSuccess+numPert.noReach;
numReachingPert = numPert.idealSuccess+numPert.eventualSuccess+numPert.noSuccess;
y_Pert_acc = [numPert.idealSuccess/numReachingPert; numPert.eventualSuccess/numReachingPert ;numPert.noSuccess/numReachingPert ;numPert.noReach/numTotalPert];
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertAccuracy.mat'],'y_Pert_acc');
if plotAccuracy == 1
    figure;
    h = bar(y_Pert_acc*100,'FaceColor','flat');
    colorScheme = colorGradient([1,0,.5],[.5,0,1],5);
    h.CData(1,:) = colorScheme(1,:);
    h.CData(2,:) = colorScheme(2,:);
    h.CData(3,:) = colorScheme(3,:);
    h.CData(4,:) = colorScheme(4,:);
    l = cell(1,4);
    l{1}='Success'; l{2}='Eventual Success'; l{3}='No Success'; l{4}='No Reach';
    set(gca,'xticklabel', l);
    ylabel('Accuracy (%)');
    title([SUB ' ' EXPER_SESSION ' ' EXPER_COND ' : Stim Accuracy']);
    %            savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajector/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
end
%% Find accuracy of eventual success /only failed reaches and multiple reaches
% Stimulation
numTotalIsolatePert = numPert.eventualSuccess + numPert.noSuccess;
y_Pert_isolate = [numPert.eventualSuccess];
y_Pert_acc_isolate = y_Pert_isolate./numTotalIsolatePert;
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_pertIsolateAccuracy.mat'],'y_Pert_acc_isolate');
% Baseline
numTotalIsolateBL = numBL.eventualSuccess + numBL.noSuccess;
y_BL_isolate = [numBL.eventualSuccess];
y_BL_acc_isolate = y_BL_isolate./numTotalIsolateBL;
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_BLIsolateAccuracy.mat'],'y_BL_acc_isolate');

%% Plot Trajectory examples (first 16 trials)
if plotSampleTrajectories == 1
    %% Plot Sample 1D Trajectories in Z
    nframes = 2500;
    % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
    plot1DTrajectories(traj1,'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials ZDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_baseline.fig']);
    plot1DTrajectories(pert.traj, 'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials ZDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_stim.fig']);
    plot1DTrajectories(ctrl.traj, 'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials ZDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_zDim_washout.fig']);
    %% Plot Sample 1D Trajectories in y
    % INPUT PARSERS AVAILABLE: coordinate_dim, startFrame, endFrame
    plot1DTrajectories(base.traj,'endFrame',nframes,'coordinate_dim',2);
    sgtitle([SUB ' ' EXPER_COND ': handpaths baseline trials YDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_baseline.fig']);
    plot1DTrajectories(pert.traj, 'endFrame',nframes,'coordinate_dim',2);
    sgtitle([SUB ' ' EXPER_COND ': handpaths stimulation trials YDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_stim.fig']);
    plot1DTrajectories(ctrl.traj, 'endFrame',nframes,'coordinate_dim',2);
    sgtitle([SUB ' ' EXPER_COND ': handpaths washout trials YDim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_yDim_washout.fig']);

    %% Plot 3D Trajectories
    plot3DTrajectories_fromTraj(traj1, 'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': 3D handpaths baseline trials']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_baseline.fig']);
    plot3DTrajectories(pert.traj, 'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': 3D handpaths stimulation trials']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_stim.fig']);
    plot3DTrajectories(ctrl.traj, 'endFrame',nframes);
    sgtitle([SUB ' ' EXPER_COND ': 3D handpaths washout trials']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SUB '_' EXPER_SESSION '_' EXPER_COND '_washout.fig']);
end
%% Find start frame of movement (uses smoothing and thresholding of velocity)
trialIdxs = [1:numTrials];
[movStrt] = getMovStartFrame(score,traj,0,trialIdxs,'allTrials',1);
%% Find end frame of reach movement (uses smoothing and thresholding of velocity)
[movEnd,~] = getMovEndFrame(score,traj,0,trialIdxs, movStrt,'allTrials',1);
%% Using JAABA - Start and End Frames of reach and grab
if useJAABA
    % Start: Hand Lift and End: Grab
    % Load classifier data (hand lift, grab, etc)
    % find how many trials there are
    TRK_SIDE =  [RAWDATA_BASEPATH SUB '/' EXPER_SESSION '/trk/'];
    foldContents = dir(TRK_SIDE);
    numTrials = 0;
    % Find how many trials there are with videos and trk files
  
    for icont = 1:numel(foldContents)
        foldEmpt = contains(foldContents(icont).name,'side');
        if foldEmpt == 0
            num2add = 0;
        else
            num2add = 1;
        end
        numTrials = numTrials + num2add;
    end
    % Initiate cells to store start frames in
    chew.t0           = {};
    atMouth.t0        = {};
    digitsTogether.t0 = {};
    lift.t0       = {};
    supinate.t0       = {};
    grab.t0           = {};
    % Initiate cells to store stop frames in
    chew.t1           = {};
    atMouth.t1        = {};
    digitsTogether.t1 = {};
    lift.t1       = {};
    supinate.t1       = {};
    grab.t1           = {};
    % load output for all trials and store in mega matrix
    for itrial = 1:numTrials
        if itrial < 10
            thisTrialNum = ['00' num2str(itrial)];
        elseif itrial >=10 && itrial < 100
            thisTrialNum = ['0' num2str(itrial)];
        elseif itrial >= 100
            thisTrialNum = num2str(itrial);
        end
        TRIAL_FOLDER = [JAABA_OUTPUT SUB '_' EXPER_SESSION '_v' thisTrialNum '/'];
        chew_scores           = load([TRIAL_FOLDER 'scores_Chewm134w.mat'],'allScores');
        atMouth_scores        = load([TRIAL_FOLDER 'scores_Atmouthm134w.mat'],'allScores');
        handOpen_scores = load([TRIAL_FOLDER 'scores_Handopenm134w.mat'],'allScores');
        lift_scores       = load([TRIAL_FOLDER 'scores_Liftm134w.mat'],'allScores');
        sup_scores       = load([TRIAL_FOLDER 'scores_Supm134w.mat'],'allScores');
        grab_scores           = load([TRIAL_FOLDER 'scores_Grabm134w.mat'],'allScores');
        % save in mega cell array start times
        chew.t0{itrial}           = chew_scores.allScores.t0s;
        atMouth.t0{itrial}        = atMouth_scores.allScores.t0s;
        digitsTogether.t0{itrial} = handOpen_scores.allScores.t0s;
        lift.t0{itrial}       = lift_scores.allScores.t0s;
        supinate.t0{itrial}       = sup_scores.allScores.t0s;
        grab.t0{itrial}           = grab_scores.allScores.t0s;
        % save in mega cell array stop times
        chew.t1{itrial}           = chew_scores.allScores.t1s;
        atMouth.t1{itrial}        = atMouth_scores.allScores.t1s;
        digitsTogether.t1{itrial} = handOpen_scores.allScores.t1s;
        lift.t1{itrial}       = lift_scores.allScores.t1s;
        supinate.t1{itrial}       = sup_scores.allScores.t1s;
        grab.t1{itrial}           = grab_scores.allScores.t1s;
    end
    %%%%JAABA OUTPUT CODE INSERT HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 % Label grab start as movEnd
    %                 grabStart.nbase = movEnd.nbase;
    %                 grabStart.npert = movEnd.npert;
    %                 grabStart.nctrl = movEnd.nctrl;

    for itrial = 1:numTrials
        if isempty(grab.t0{1,itrial}{1,1})
            grabStart = 0;
        else
            grabStart(1,itrial) = grab.t0{1,itrial}{1,1}(1);
        end
    end

end % JAABA
%% Get num of success, no success, etc for each baseline,stim, and wash
% Score each trial
scores = NaN(1,1:numTrials);
for itrial = 1:numTrials
    this_chew = chew.t0{1,itrial}{1,1}(1);
    this_lift = lift.t0{1,itrial}{1,1}(1);
    scores(1,itrial) = 
end
% baseline
numBL = getNumTrialsPerScore(trialIdxs, base);
% stim
numPert = getNumTrialsPerScore(trialIdxs, pert);
% washout
numWash = getNumTrialsPerScore(trialIdxs, ctrl);
%% Plot Overlay 1D - only reach
if plot1DOverlayTrajectories == 1
    % Define deminsion to plot (x,y,z)
    dimPlot = 3; % if change add dimPlot input to plotOverlay1D function
    % Plot overlay for baseline
    movIdx = [1:numTrials];
    plotOverlay1D(movStrt, movEnd,movIdx,traj);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline_dim' num2str(dimPlot) '.fig']);
    % Plot overlay for stimulation
    plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj);
    title([SUB ' ' EXPER_COND ': Overlay' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachStim_dim' num2str(dimPlot) '.fig']);
    % Plot overlay for all baseline, pre and post stimulation
    movStrtBL = [movStrt.nbase movStrt.nctrl];
    movEndBL = [movEnd.nbase movEnd.nctrl];
    movStrtIdxs = [movIdx.nbase movIdx.nctrl]';
    plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE 'pre and post stim (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachAllBL_dim' num2str(dimPlot) '.fig']);
    %% Plot Overlay 1D - whole movement
    % Plot baseline
    plotOverlay1D(movStrt, movEnd,movIdx, traj,'PlotEntireReach',1);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayBaseline_dim' num2str(dimPlot) '.fig']);
    % Plot stimulation
    plotOverlay1D(movStrt.npert, movEnd.npert,movIdx.npert,traj,'PlotEntireReach',1);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' stim Trials (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayStim_dim' num2str(dimPlot) '.fig']);
    % Plot all baseline, pre and post stimulation
    movStrtBL = [movStrt.nbase movStrt.nctrl];
    movEndBL = [movEnd.nbase movEnd.nctrl];
    movStrtIdxs = [movIdx.nbase movIdx.nctrl]';
    plotOverlay1D(movStrtBL, movEndBL,movStrtIdxs,traj,'PlotEntireReach',1);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim (Dim ' num2str(dimPlot) ')']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/1DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayAllBL_dim' num2str(dimPlot) '.fig']);
end
%% Plot Overlay 3D - reach movement
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_Overlay3DVariables.mat'],'movStrt', 'movEnd','movIdx', 'traj');
if plot3DOverlayTrajectories == 1
    % Plot baseline
    plotOverlay3D(movStrt, movEnd,movIdx, traj);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' baseline Trials']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE '_overlayReachBaseline.fig']);
    % Plot all baseline, pre and post stimulation
    movStrtBL = [movStrt.nbase movStrt.nctrl];
    movEndBL = [movEnd.nbase movEnd.nctrl];
    movStrtIdxs = [movIdx.nbase movIdx.nctrl]';
    plotOverlay3D(movStrtBL, movEndBL,movStrtIdxs,traj);
    title([SUB ' ' EXPER_COND ': Overlay ' SCORE ' pre and post stim']);
    savefig([BASEPATH 'Figures/' SUB '/Behavior/3DTrajectories/' SCORE '/' SUB '_' EXPER_SESSION '_' EXPER_COND '_' SCORE 'overlayReachAllBL.fig']);
end
%% Get path length for baseline, stim, and washout trials
% Get path lengths for baseline trials
[totalPathLength] = getPathLength(traj, movIdx,movStrt, movEnd);
% Get path lengths for stimulation trials
[totalPathLength.npert] = getPathLength(traj, movIdx.npert,movStrt.npert, movEnd.npert);
% Get path lengths for washout baseline trials
[totalPathLength.nctrl] = getPathLength(traj, movIdx.nctrl,movStrt.nctrl, movEnd.nctrl);
% Save path lengths in mat file
save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_PathLengths.mat'],'totalPathLength');
%% Get digit spread
if useJAABA
    % get pointer digit at first grab for each trial
    % traj(trials x coordinate xyz x frames)
    % find traj points for start of grab
    trajDig1 = traj(movIdx,:,grabStart);
    trajDig2 = traj2(movIdx,:,grabStart);

    % baseline trial spread
    for i = 1:length(movIdx)
        dig1 = trajDig1(i,:,i);
        dig2 = trajDig2(i,:,i);
        digitSpread(i) = sqrt(sum((dig1-dig2) .^ 2));
    end
    save([BASEPATH 'Data_Analyzed/' SUB '/Behavior/' SUB '_' EXPER_SESSION '_' EXPER_COND '_DigitSpread.mat'],'digitSpread','movIdx','trialIdxs');
end

%% Compare accuracy across conditions
Behavioral_Metrics(animals,BASEPATH,exper_conditions);
%% Compare 3D trajectories
%Behavioral_Comparison_Conditions; %(animals,BASEPATH, RAWDATA_BASEPATH,exper_conditions,SCORE,USER);
%% Compare digit spread (pointer and ring), harm vs control
% Calculate distance between digits
% at grab start using traj
Behavioral_DigitSpread; %trajDigSpread, movIdx, trialIdxs);
%% Visualize liftTime
