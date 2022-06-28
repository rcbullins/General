function create_trial_structure_rcb(rootdir,neural_dir,wasabi_vid,save_dir,...
    recsite,expt_type)

plotFig=0;
saveData=1;
%
% neural_data_dir = neural_data_dirs{i}
% vid_data_dir = vid_data_dirs{i}
% expt_type = exp_type(i)
% recsite = site{i}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the timestamps and create session trial labels
clear con_trial LT_only LT_single LT_4hz LT_10hz LT_20hz LT_40hz
clear sync table_start trial_start laser laser_start tone light
clear laser_gate laser_gate_start trial_array

% load([neural_data_dir '/event_ind.mat']);
load([neural_dir filesep 'event_ind.mat']);

% get recording sample rate
% mname = ls([neural_data_dir filesep recsite '\*.nidq.meta']);
mname = ls([neural_dir filesep recsite '\*.nidq.meta']);
% metafile = [neural_data_dir filesep recsite filesep mname];
metafile = [neural_dir filesep recsite filesep mname(end,:)];
meta_text = fileread([metafile]);
[foo bar] = regexp(meta_text,'niSampRate=\d');
mystr = meta_text(bar+(0:5));
if ~isempty(str2num(mystr(end)))
    fs = str2num(mystr);
else
    fs = str2num(mystr(1:(end-1)));
end

if plotFig
    close all;
    figure();
    set(gcf,'Position',[1850 20 1500 1300]);
    
    subplot(6,1,1)
    hold on
    %         plot(sync,2e-3*ones(length(sync),1)','k.')
    title('sync')
    
    subplot(6,1,2)
    hold on
    plot(trial_start,2e-3*ones(size(trial_start,2),1),'g*')
    title('camera')
    
    subplot(6,1,3)
    hold on
    plot(table_start,2e-3*ones(size(table_start,2)),'r*')
    title('table_start')
    
    subplot(6,1,4)
    plot(laser,2e-3*ones(size(laser,2)),'c*')
    hold on
    % laser pulse timestamps
    plot(laser_start,2e-3*ones(size(laser_start,2)),'m*')
    title('laser')
    
    if expt_type==1
        subplot(6,1,5)
        hold on
        plot(tone,2e-3*ones(size(tone,2)),'r*')
        title('tone')
        
        if ~isempty(light)
            subplot(6,1,6)
            hold on
            plot(light,-2e-3*ones(size(light,2)),'b*')
            title('masking light')
        end
    elseif expt_type==2
        
        % laser pulse timestamps
        subplot(6,1,5)
        plot(laser_gate,2e-3*ones(size(laser,2)),'c*')
        hold on
        plot(laser_gate_start,2e-3*ones(size(laser_start,2)),'m*')
        title('laser_gate')
    end
    
    zoom on;
    
    %     linkaxesInFigure('x');
end

frontList = dir(fullfile([wasabi_vid '\tracked\'], '*front*.trk'));
sideList = dir(fullfile([wasabi_vid '\tracked\'], '*side*.trk'));
disp(['front_vid: ' num2str(numel(frontList)) ' side_vid: ' num2str(numel(sideList)) ', trial_start: ' num2str(length(trial_start))])

if numel(frontList) ~= numel(sideList) | numel(frontList) ~= length(trial_start)
    disp(['!!! Number of trial_start is NOT equal to number of videos !!! '])
else
    disp(['Number of trial_start is equal to number of videos :) '])
end
if strcmp(neural_dir,[rootdir '\jcr79\ephys\jcr79_20200930_4500_915um_2_g0'])
    vid_offset=28;
    disp(['front_vid: ' num2str(numel(frontList)-vid_offset) ' side_vid: ' num2str(numel(sideList)-vid_offset) ', trial_start: ' num2str(length(trial_start))])
end
if strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0'])
    vid_offset=-3;
    disp(['front_vid: ' num2str(numel(frontList)-vid_offset) ' side_vid: ' num2str(numel(sideList)-vid_offset) ', trial_start: ' num2str(length(trial_start))])
end
disp([neural_dir]);
disp([wasabi_vid]);


% % create 6 column array for all trials
% (t,1) Use trial flag
% (t,2) trial start time
% (t,3) laser trial flag
% (t,4) number of pulses in stimulus train
% (t,5) frequency of stimulus train
% (t,6) laser_gate trial flag
trial_array=[]; % clear trial array
% description of list of variables in trial array
trial_array_var={'use/exclude trial flag' 'table cue start' 'laser start' 'number of laser pulses' 'laser frequency' 'laser_gate_start'};

for t = 1:length(trial_start)
    if expt_type==1
        if strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0']) && (t==103 | t==104 | t==105)
            trial_array(t,1) = 0; % 0/1 flag to use trial
        else
            trial_array(t,1) = 1; % 0/1 flag to use trial
        end
        tab = table_start(table_start>=trial_start(t)-5 & table_start<=trial_start(t)+5e4); % flag for table_start/cue
        if tab
            trial_array(t,2) = tab;
        end
        
        las = laser(laser>=trial_start(t)-5 & laser<=(trial_start(t)+10e4)); % flag for laser trial
        if las
            trial_array(t,3) = las(1);
        end
        
        trial_array(t,4) = length(las); % numPulses in laser train
        if trial_array(t,4)>1
            trial_array(t,5) = round(1e6/mean(diff((1e6/fs)*las))); % in Hz
        end
        
    elseif expt_type==2
        
        trial_array(t,1) = 1; % 0/1 flag to use trial
        
        tab = table_start(table_start>=trial_start(t)-5 & table_start<=trial_start(t)+5e4); % flag for table_start/cue
        if tab
            trial_array(t,2) = tab;
        end
        
        las = laser(laser>=trial_start(t)-5 & laser<=(trial_start(t)+10e4)); % flag for laser trial
        if las
            trial_array(t,3) = las(1);
        end
        
        trial_array(t,4) = length(las); % numPulses in laser train
        if trial_array(t,4)>1
            trial_array(t,5) = round(1e6/mean(diff((1e6/fs)*las))); % in Hz
        end
        
        las_gate = laser_gate(laser_gate>=trial_start(t)-5 & laser_gate<=(trial_start(t)+10e4)); % flag for laser trial
        if las_gate
            trial_array(t,6) = las_gate(1);
        end
        elseif expt_type==3
        
        trial_array(t,1) = 1; % 0/1 flag to use trial
        
        tab = table_start(table_start>=trial_start(t)-5 & table_start<=trial_start(t)+5e4); % flag for table_start/cue
        if tab
            trial_array(t,2) = tab;
        end
        
        las = laser(laser>=trial_start(t)-5 & laser<=(trial_start(t)+10e4)); % flag for laser trial
        if las
            trial_array(t,3) = las(1);
        end
        
        trial_array(t,4) = length(las); % numPulses in laser train
        if trial_array(t,4)>1
            trial_array(t,5) = round(1e6/mean(diff((1e6/fs)*las))); % in Hz
        end
        
        
    end
end

% categorize trials, create lists of trial_type indices
if expt_type==1
    con_trial = find(trial_array(:,5)==0 & trial_array(:,2)>0); % con_trial trial_start
    LT_single = find(trial_array(:,2)>0 & trial_array(:,4)==1); % light single pulse
    LT_only = find(trial_array(:,5)==10 & trial_array(:,2)==0); % light only, 10hz
    LT_4hz = find(trial_array(:,5)==4 & trial_array(:,2)>1); % cue and light, 4hz
    LT_10hz = find(trial_array(:,5)==10 & trial_array(:,2)>1); % cue and light, 10hz
    LT_40hz = find(trial_array(:,5)==40 & trial_array(:,2)>1); % cue and light, 40hz
    LT_gate=[];
elseif expt_type==2
    % for exp_type==2, jay's data
    % There are 3 laser trial types to categorize
    % NOTE: the laser_gate channel has signal for both trial types
    % 3 trial types:
    % 1== if CUE && LASER && LASER_GATE onset are within +/- 5 timestamps, 
    % then its a 'nHz laser trial'
    % 2== if CUE && LASER onset are within +/- 5 timestamps && LASER_GATE onset > CUE & LASER, 
    % then its a 'nHz laser_gate trial'
    % 3== if CUE && LASER onset are within +/- 5 timestamps && no LASER_GATE (==0), 
    % then its a 'missed LASER_GATE trial' and is a 'control trial'
    
    con_trial = find(((trial_array(:,3)==0 & trial_array(:,6)==0) | ...
        (trial_array(:,6)>0 & trial_array(:,3)==0) | ...
        (trial_array(:,6)==0 & trial_array(:,3)>0)) & ...
        trial_array(:,2)>0); % con_trial trial_start
    LT_only = find(trial_array(:,5)==10 & trial_array(:,6)>0 & trial_array(:,2)==0); % laser only
    LT_4hz = find(trial_array(:,5)==4 & trial_array(:,6)>0 & trial_array(:,2)>1); % cue and light, 4hz
    LT_10hz = find(trial_array(:,5)==10 & trial_array(:,6)>0 & trial_array(:,2)>1); % cue and light, 10hz
    LT_40hz = find(trial_array(:,5)==40 & trial_array(:,6)>0 & trial_array(:,2)>1); % cue and light, 40hz
    LT_gate = find((trial_array(:,6)-trial_array(:,3))>-5); % laser_gate, 40hz
elseif expt_type==3
    % for exp_type==3, jay's data, whisper, no sync pulses
    % There are 2 laser trial types to categorize, laser/cue and laser only
    % 2 trial types:
    % 1== if CUE && LASER onset are within +/- 5 timestamps,
    % then its a 'nHz laser trial'
    
    con_trial = find(trial_array(:,5)==0 & trial_array(:,2)>0); % con_trial trial_start
    LT_only = find(trial_array(:,5)==10 & trial_array(:,2)==0); % light only, 10hz
    LT_4hz = find(trial_array(:,5)==4 & trial_array(:,2)>1); % cue and light, 4hz
    LT_10hz = find(trial_array(:,5)==10 & trial_array(:,2)>1); % cue and light, 10hz
    LT_40hz = find(trial_array(:,5)==40 & trial_array(:,2)>1); % cue and light, 40hz
    LT_gate=[];
    
end

if saveData
%     save([save_dir filesep 'event_ind.mat'],'trial_array','con_trial','LT_only',...
%         'LT_4hz','LT_10hz','LT_40hz','LT_gate','-append');
     
    save([save_dir filesep 'event_ind.mat'],'trial_array','con_trial',...
        'LT_only','LT_single','LT_4hz','LT_10hz','LT_40hz','LT_gate','-append');
end

disp([neural_dir ' - create trial structure complete'])
