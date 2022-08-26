% Step 3: Pre-processing Neural and Behavior Data in to Struct.
%
% Processing of Motor Skill Reaching Task.
% 1) extract pulse timestamps from recording file (eg .nidq or .imec)
%%% camera triggers, lights, lasers, tables, tones, etc...
% 2) create variables that describe the experiment trial structure
%%% control trials, light trials, light only trials, etc...
% 3) create variables that describe trial by trial outcomes
%%% use/exclude trial flag, trial time, success/fail/, lift, reach,
%%% reaction time, reach velocity,

% Previoulsy: jcr_proc_log.m or rcb_proc_log
%% Setup Project Names
ProjectName = 'ThalamusSP_Project';
SUB = 'jcr70';
EXPER_SESSION = '20191016';

%% Add Pathways
BASEPATH = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\';
PROJECT  = [BASEPATH ProjectName '\'];
CALIB = [PROJECT 'Data\Calib_Results_stereo.mat'];

CODE_SP = [PROJECT 'Code\'];
CODE_GENERAL = [BASEPATH 'Code\General\'];
CODE_CALIB = [PROJECT 'Data\Calib_Results_stereo.mat'];
CODE_BRITTON = [BASEPATH 'Code\Packages\britton_code\code\'];
CODE_BRITTON_PLOT = [BASEPATH 'Code\Packages\britton_code\other_code'];
CODE_PROBES = [BASEPATH 'Code\Packages\Probes\'];
SPIKE_GLX = [BASEPATH 'Code\Packages\SpikeGLX_Datafile_Tools\'];

addpath(genpath(CODE_GENERAL));
addpath(CODE_CALIB);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));
addpath(genpath(CODE_PROBES));
addpath(genpath(SPIKE_GLX));
%% Set Directories of Experimental List
% Video directories
VID_DATA = {[PROJECT 'Data\jcr70\video\20191016']};

% Corresponding .nidq event ind directories
NEURAL_DATA = {[PROJECT 'Data\jcr70\ephys\jcr70_20191016_920um_g0']};

% Corresponding JAABA Dir
JAABA_OUTPUT = {[PROJECT 'Data_Analyzed\jcr70\JAABA_Output\jcr70_20191016_JAABA_Output.mat']};

% subfolder location of the nidq meta file, which states the number of channels
site = {'cortex'};

% experiment type (eg, jcr or jay)
% 1==jcr jeremy rig, 2==jay janelia new rig, 3==jay janelia old rig
exp_type = [1];

dir_list = 1; %[1:2 5 7 8:numel(neural_data_dirs)]; % list of files to process
%% Define Channels in nidq
% JCR Experiments, Jeremy Cohen
% timestamps on Jeremy's janelia rig in .nidq file
% Get data from .nidq file channel #
% for 'exp_type==1'
ch1=64;% sync pulse timestamps, jc_pulse 1hz
ch2=66;% camera timestamps
ch3=68;% table_start/cue pulse timestamps
ch4=67;% laser pulse timestamps
ch5=70;% tone pulse timestamps
ch6=69;% masking led light pulse timestamps
jcr_pulse_channels = [ch1 ch2 ch3 ch4 ch5 ch6];

% Jay Guo Experiments
% timestamps on Jay's janelia corner rig in .nidq file
% Get data from .nidq file channel #
% for 'exp_type==2'
ch1=0;% sync pulse timestamps, imec 1hz
ch2=1;% camera timestamps
ch3=2;% table_start/cue pulse timestamps
ch4=3;% laser pulse timestamps
ch5=4;% laser gate, for grab detector trials
jay_pulse_channels_1 = [ch1 ch2 ch3 ch4 ch5];

% for 'exp_type==3'
ch1=64;% camera timestamps
ch2=65;% table_start/cue pulse timestamps
ch3=67;% laser pulse timestamps
ch4=66;% laser gate, for grab detector trials
jay_pulse_channels_2 = [ch1 ch2 ch3 ch4];
%% Loop through each experiment in list
smoothwin = 25; % smoothing window (for plotting purposes)
tic;
% For each session
for i=[1]%dir_list(5:end)
    % set rootdir, neural data, and video data directories
        this_NEURAL = NEURAL_DATA{i};
        this_VID = VID_DATA{i};
        this_TRK = [this_VID 'tracked\']; 
        this_JAABA = JAABA_OUTPUT{i};
        this_SAVE = this_NEURAL;
    
    % set channels for event processing
    if exp_type(i)==1
        channels=jcr_pulse_channels;
    elseif exp_type(i)==2
        channels=jay_pulse_channels_1;
    elseif exp_type(i)==3
        channels=jay_pulse_channels_2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) Extract the timestamps for every channel in .nidq file
    % stored in the local directory in a file '/ind.mat'
   
     get_event_ind_rcb(this_NEURAL,site{i},channels,exp_type(i),this_SAVE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) Set up trial array and categorize trials into trial type
    % append new variables to the event ind.mat file
     create_trial_structure_rcb(PROJECT,this_NEURAL,this_TRK,this_SAVE,...
         site{i},exp_type(i))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) Create trial outcome variables, lift, reach, grasp, handopen, 
    % supinate, at mouth, success or fail outcome, trajectories.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Experiment specific info - DO NOT CHANGE %
    % TO ALIGN LASER+CUE AND CONTROL TRIALS, align to pellet shift onset
    if exp_type(i)==1
        if i < 9
            % FOR JCR80_20201103 and earlier experiments
            frames = 1:2000; % frame range (2ms frame rate)
            LT_shift = 250; % frame shift LT trials
            CON_shift = 25;
        else
            % FOR JCR80_20201104 and later experiments
            frames = 1:2000; % frame range (2ms frame rate)
            LT_shift = 25; % frame shift LT trials
            CON_shift = 25;
        end
    elseif exp_type(i)==2 || exp_type(i)==3
        frames = 1:2000; % frame range (2ms frame rate)
        LT_shift = 0; % frame shift LT trials
        CON_shift = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotfig=0; % plot trial trajectories during processing
    create_trial_outcome_rcb(PROJECT,this_NEURAL,this_VID,this_SAVE,...
         CALIB,exp_type(i),smoothwin,frames,LT_shift,CON_shift,plotfig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) Get LFP chunks and add to struct
    fname = ls([this_NEURAL filesep site{i} filesep '*.nidq.bin*']);
    RAW_EPHYS_FILE = ([this_NEURAL filesep site{i} filesep fname]);

    lfp = getLFPfromBin(RAW_EPHYS_FILE,BASEPATH, SUB, EXPER_SESSION);
    saveData = 1;
    if saveData
        % note some variables not listed here may be previously saved in
        % 'create_trial_structure_jcr.m'
        save([this_SAVE '\lfp.mat'], ...
          'lfp','-v7.3');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5) create the final data structure with all processed data
    % neural data and kinematics combined.
    new_or_append = 0; %==1 update/append to current file, ==0 create new file
    if new_or_append==0
        event_ind_dir = NEURAL_DATA{i};
        ctx=1; % 1==process data, 0==do not process
        thal=0; % 1==process data, 0==do not process
        create_processed_data_struct_rcb(event_ind_dir,JAABA_dir,ctx,thal)
    elseif new_or_append==1
        update_processed_data_struct(event_ind_dir)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
b=toc;
disp(b);

%% PLOT - Individual Trial Position Data
t_win=1:4000; % time window in ms
smoothwin = 25; % 25 frames, 50ms
event_dir = [NEURAL_DATA{i} filesep 'event_ind']; % event_ind.m location
trial=1:5;
plot_trial_rcb(event_dir,trial,t_win,smoothwin,exp_type(i));
