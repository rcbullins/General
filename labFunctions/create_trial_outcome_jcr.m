function create_trial_outcome_jcr(rootdir,neural_dir,vid_dir,save_dir,...
    calib_file,expt_type,smoothwin,frames,LT_shift,CON_shift,plotfig)

% create trial outcome variables based on derived behavioral timepoints
% first, load 3D hand trajectories based on tracked video files (APT) and
% the stereo_calibration.mat file.
% then create behavioral timepoints: lift, reach, grasp, atmouth, etc...
% then create trial outcome variables, such as: success_rate, reach_rate,
% reaction_time, etc...

% expt_type=exp_type(i);
%%
saveData=1;

disp(neural_dir)
disp(vid_dir)
disp(save_dir)
nframes=frames(end); % frames to grab
tic;
if expt_type==1
    i_obj = 1; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_1 conf_1] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 2; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_2 conf_2] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 3; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_3 conf_3] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 4; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_4 conf_4] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
elseif expt_type==2
    i_obj = 1; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_1 conf_1] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 2; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_2 conf_2] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 3; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_3 conf_3] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 4; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_4 conf_4] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
elseif expt_type==3
    i_obj = 1; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_1 conf_1] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 2; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_2 conf_2] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 3; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_3 conf_3] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
    i_obj = 4; % object number to track (1=2nd digit, 2=4th digit, 3=wrist, 4=pellet)
    [traj_4 conf_4] = get_traj_3D([vid_dir filesep 'tracked' filesep],calib_file,nframes,i_obj);
end


toc;
% clear variables
clear con_trial LT_only LT_4hz LT_10hz LT_40hz
clear sync table_start tone trial_start laser laser_start light trial

load([neural_dir filesep 'event_ind.mat']);

% to account for split files or missed video
vid_offset=0;
if strcmp(neural_dir,[rootdir '\jcr79\ephys\jcr79_20200930_4500_915um_2_g0'])
    vid_offset=28;
    frontList = dir(fullfile([vid_dir '\tracked\'], '*front*.trk'));
    sideList = dir(fullfile([vid_dir '\tracked\'], '*side*.trk'));
    disp(['front_vid: ' num2str(numel(frontList)-vid_offset) ' side_vid: ' num2str(numel(sideList)-vid_offset) ', trial_start: ' num2str(length(trial_start))])
    
elseif strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0'])
    vid_offset=-3;
    frontList = dir(fullfile([vid_dir '\tracked\'], '*front*.trk'));
    sideList = dir(fullfile([vid_dir '\tracked\'], '*side*.trk'));
    disp(['front_vid: ' num2str(numel(frontList)) ' side_vid: ' num2str(numel(sideList)) ', trial_start: ' num2str(length(trial_start))])
    disp(['Trial 102 video contains 10010 frames for four trial_start, 102:105'])
end


% create and clear behavioral variables
% Lift times and trial outcome
digvel=[];digaccel=[];wvel=[];waccel=[];accelvec=[];
digdist3d=[]; d_digdist3d=[];
digpelletdist3d=[]; trial_outcome=[]; reaction_time=[];
dist1d_x =[]; dist1d_y=[]; dist1d_z=[]; dist3d=[];
reach_vel=[];
lift_flag=[];
lift_t=nan(length(trial_start),frames(end));
lift_smoothwin=smoothwin;

if strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0'])
    trial_range=[1:102 106:length(trial_start)];
elseif strcmp(neural_dir,[rootdir '\M341\ephys\M341_20210817_1000_g0_imec0'])
    trial_range=[1:120];
else
    trial_range=1:length(trial_start);
end

trial_outcome_array = [];
% trial_outcome_array(t,:) = [trial_outcome(t) reach_flag(t) lift_flag(t) reach_vel(t) reaction_time(t) lift(1) reach(1) grasp(1) pelletgrasp(1) digmouth(1) pelletmouth(1)];
% 1) trial_outcome; % success==1, fail==0
% 2) lift_flag; % lift==1, no lift==0
% 3) reach_flag; % reach==1, no reach==0
% 4) reach_vel; % distance from lift to grab / time
% 5) reaction_timee; % table/cue start to lift time
% 6) lift; % success==1, fail==0
% 7) reach; % success==1, fail==0
% 8) grasp; % success==1, fail==0
% 9) pelletgrasp; % success==1, fail==0
% 10) digmouth; % success==1, fail==0
% 11) pelletmouth; % success==1, fail==0

% acceleration threshold
if expt_type==1
    if strcmp(neural_dir,'D:\localData\jcr70\ephys\jcr70_20191016_920um_g0')
        accel_thresh=0.9e4; % 
    else
        accel_thresh=1e4; % 
    end
elseif expt_type==2 || expt_type==3
    accel_thresh=1.2e4; % 
end
% the crossing of this value is used to set reach and lift indices.

tic;

% plotfig=0;
for t = trial_range
    disp(['trial ' num2str(t)])
    
    % clear these behavioral variables for each trial
    clear lift reach grasp pelletgrasp digmouth pelletmouth
    
    if strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0']) && t>102
        vid_offset=-3;
    else
        vid_offset=0;
    end
    
    % dig2
    d2x = smoothdata(squeeze(traj_1(t+vid_offset,1,:))','movmean',lift_smoothwin);
    d2y = smoothdata(squeeze(traj_1(t+vid_offset,2,:))','movmean',lift_smoothwin);
    d2z = smoothdata(squeeze(traj_1(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % dig4
    d4x = smoothdata(squeeze(traj_2(t+vid_offset,1,:))','movmean',lift_smoothwin);
    d4y = smoothdata(squeeze(traj_2(t+vid_offset,2,:))','movmean',lift_smoothwin);
    d4z = smoothdata(squeeze(traj_2(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % wrist
    wx = smoothdata(squeeze(traj_3(t+vid_offset,1,:))','movmean',lift_smoothwin);
    wy = smoothdata(squeeze(traj_3(t+vid_offset,2,:))','movmean',lift_smoothwin);
    wz = smoothdata(squeeze(traj_3(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % pellet
    px = smoothdata(squeeze(traj_4(t+vid_offset,1,:))','movmean',lift_smoothwin);
    py = smoothdata(squeeze(traj_4(t+vid_offset,2,:))','movmean',lift_smoothwin);
    pz = smoothdata(squeeze(traj_4(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    
    % vel, accel, and distances
    for f=frames
        if f>1
            % digit movement
            digdist3d(t,f) = norm([d4x(f) d4y(f) d4z(f)] - [d2x(f-1) d2y(f-1) d2z(f-1)]); % inter-digit distance
            d_digdist3d(t,f) = digdist3d(t,f) - digdist3d(t,f-1); % change in digit distance, to capture grasp onset
            digvel(t,f) = (d4x(f)-d4x(f-1))/0.002; % digit-4 velocity
            digaccel(t,f) = (digvel(t,f) - digvel(t,f-1))/0.002; % digit-4 accel
            digpelletdist3d(t,f) = norm([d4x(f) d4y(f) d4z(f)] - [px(f) py(f) pz(f)]);
            
            % wrist vel acceleration
            wvel(t,f) = (wx(f)-wx(f-1))/0.002;
            waccel(t,f) = (wvel(t,f)-wvel(t,f-1))/0.002;
        end
    end
    
    % create accel vector, composed of dig + wrist accel
    accelvec(t,:) = (abs(digaccel(t,:)))+(abs(waccel(t,:)));
    accelvec_sm(t,:) = smoothdata(accelvec(t,:),'movmean',5);
    % find  indices that cross threshold
    accel_ind = find(accelvec(t,1:(end-1))<accel_thresh & accelvec(t,2:end)>=accel_thresh);
    accel_sm_ind = find(accelvec_sm(t,1:(end-1))<=accel_thresh & accelvec_sm(t,2:end)>=accel_thresh);
    % set start trial frame rate and range. start to +1000 frames
    if expt_type==1
        if LT_shift==250 && CON_shift==25 && ismember(t,con_trial)
            outcome_frames=500:1750; % con trial start frame for all con_trial trial_start
        elseif LT_shift==25 && CON_shift==25 && ismember(t,con_trial)
            outcome_frames=250:1500; % con trial start frame for all con_trial trial_start
        elseif strcmp(neural_dir,'D:\localData\jcr70\ephys\jcr70_20191016_920um_g0')
            outcome_frames=250:2000;
        end
    elseif expt_type==2 || expt_type==3
        outcome_frames=500:1500; % LT trial start frame for these experiments
    end
    
    % acceleration vector (accelvec) plot, and position.
    if plotfig
        close all;
        figure; set(gcf,'position',[1800 700 1500 600],'Name',neural_dir);
        hold on;
        sp(1)=subplot(7,1,1); plot(d4x); ylabel('d4x');
        sp(2)=subplot(7,1,2); plot(d4y); ylabel('d4y');
        sp(3)=subplot(7,1,3); plot(d4z); ylabel('d4z');
        sp(4)=subplot(7,1,4); plot(wx); ylabel('wx');
        sp(5)=subplot(7,1,5); plot(wz); ylabel('wz');
        sp(6)=subplot(7,1,6); plot(pz); ylabel('pz');
        sp(7)=subplot(7,1,7);
        plot(accelvec(t,:)); hold on; plot(xlim,[accel_thresh accel_thresh],'k--');
        plot(accel_ind,accel_thresh*ones(length(accel_ind),1)','m*')
        plot(accelvec_sm(t,:)); hold on; plot(xlim,[accel_thresh accel_thresh],'--','color',[0.7 0.7 0.7]);
        plot(accel_sm_ind,accelvec_sm(t,accel_sm_ind),'r*')
        
        set(gca,'xlim',[outcome_frames(1) outcome_frames(end)])
        ylabel('accel')
        linkaxes(sp,'x');
        zoom on
    end
    
    % outcome
    if expt_type==1
        lift = find(accelvec_sm(t,outcome_frames)>(accel_thresh));
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.9));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.8));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.7));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.6));
        end
        
        if isempty(lift)
            lift=nan;
        else
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+lift(1),d4x(1,outcome_frames(1)+lift(1)),'L')
                subplot(7,1,3)
                text(outcome_frames(1)+lift(1),d4z(1,outcome_frames(1)+lift(1)),'L')
                subplot(7,1,5)
                text(outcome_frames(1)+lift(1),wz(1,outcome_frames(1)+lift(1)),'L')
                text(outcome_frames(1),max(wz)*0.9,['Lift = ' num2str(outcome_frames(1)+ lift(1))])
            end
        end
        
        if strcmp(neural_dir,'D:\localData\jcr70\ephys\jcr70_20191016_920um_g0')
            reach = find(d4x(outcome_frames)>6 & smoothdata(accelvec(t,outcome_frames),'movmean',smoothwin)>accel_thresh/2);    
        else
            reach = find(d4x(outcome_frames)>-5 & smoothdata(accelvec(t,outcome_frames),'movmean',smoothwin)>accel_thresh/2);
        end
        if isempty(reach)
            reach=nan;
        else
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+reach(1),d4x(1,outcome_frames(1)+reach(1)),'R')
                subplot(7,1,3)
                text(outcome_frames(1)+reach(1),d4z(1,outcome_frames(1)+reach(1)),'R')
                subplot(7,1,5)
                text(outcome_frames(1)+reach(1),wz(1,outcome_frames(1)+reach(1)),'R')
            end
        end
        
        % grasp
        dist3d = d_digdist3d(t,:);
        dmax = find(dist3d==max(dist3d(outcome_frames(1):outcome_frames(1)+300)));
        dmin = find(dist3d==min(dist3d(dmax:outcome_frames(end))));
        grasp_f_all_ind = find(dist3d(dmax:dmin)*100<-1);
        grasp_f = grasp_f_all_ind + dmax;
        if isempty(grasp_f)
            grasp=nan;
        else
            grasp = grasp_f(1)-2;
            if plotfig
                subplot(7,1,3)
                text(grasp(1),d4z(1,grasp(1)),'G')
                subplot(7,1,6)
                text(grasp(1),pz(1,grasp(1)),'G')
            end
        end
%         if plotfig
%             figure()
%             hold on;
%             plot(outcome_frames,dist3d(1,outcome_frames).*100);
%             if lift
%                 
%                 plot(x,dist3d(1,x).*100);
%             end
%             title('Inter-digit Distance and Change in Distance')
%             xlabel('time (ms)')
%         end
        pelletgrasp = find(digpelletdist3d(t,outcome_frames)<=3);
        if isempty(pelletgrasp)
            pelletgrasp=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+pelletgrasp(1),pz(1,outcome_frames(1)+pelletgrasp(1)),'P')
                subplot(7,1,4)
                text(outcome_frames(1)+pelletgrasp(1),pz(1,outcome_frames(1)+pelletgrasp(1)),'P')
            end
        end
        
        digmouth = find((d4x(outcome_frames)>5 & d4x(outcome_frames)<15) & d4y(outcome_frames)>170 & d4z(outcome_frames)>-1);
        if isempty(digmouth)
            digmouth=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+digmouth(1),d4z(1,outcome_frames(1)+digmouth(1)),'DM')
            end
        end
        pelletmouth = find((px(outcome_frames)>5 & px(outcome_frames)<15) & py(outcome_frames)>170 & pz(outcome_frames)>0);
        if isempty(pelletmouth)
            pelletmouth=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+pelletmouth(1),d4z(1,outcome_frames(1)+pelletmouth(1)),'PM')
                subplot(7,1,4)
                text(outcome_frames(1)+pelletmouth(1),pz(1,outcome_frames(1)+pelletmouth(1)),'PM')
            end
        end
        
        % hand open
        % interdigit distance
        dist3d = d_digdist3d(t,:);
        hopen_f = find(dist3d(outcome_frames)*100>4);
        if ~isempty(hopen_f)
            hopen = hopen_f(1)+outcome_frames;
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+hopen(1),d4x(1,outcome_frames(1)+hopen(1)),'H')
            end
        else
            hopen = nan;
        end
        
        
        %supination
        % hand angle
        slope = (d2y - d4y) ./ (d2z - d4z);
        angle = atand(slope);
        angle(angle>0)=90-angle(angle>0);
        angle(angle<0)=90-angle(angle<0);
        ang = smoothdata(angle,'movmean',smoothwin);
        d_ang=[0];
        for a=1:length(ang)
            if a>1
                d_ang(a) = ang(a) - ang(a-1);
            else
                d_ang(a)=0;
            end
        end
        if ~isnan(grasp)
        dmin = find(d_ang==min(d_ang(grasp(1)-10:grasp(1)+50)));
        end
        if ~isempty(dmin)
            supinate = dmin;
            if plotfig
                subplot(7,1,5)
                text(supinate(1),wz(1,supinate),'S')
            end
        else
            supinate = nan;
        end
                
        if plotfig
            figure; set(gcf,'position',[1800 110 1500 500]);
            hold on;
            plot(outcome_frames,dist3d(outcome_frames).*100);
            title('Change in Inter-digit Distance,grasp=star')
            xlabel('frame')
            if grasp>0
                plot(grasp(1),dist3d(grasp(1))*100,'m*')
                text(outcome_frames(1)+50,max(dist3d(outcome_frames))*100*0.9,['grasp = ' num2str(grasp(1))]);
            end
        end
    elseif expt_type==2 || expt_type==3
        lift = find(accelvec_sm(t,outcome_frames)>(accel_thresh));
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.9));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.8));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.7));
        end
        if isempty(lift)
            lift = find(accelvec(t,outcome_frames)>(accel_thresh*0.6));
        end
        
        if isempty(lift)
            lift=nan;
        else
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+lift(1),d4x(1,outcome_frames(1)+lift(1)),'L')
                subplot(7,1,3)
                text(outcome_frames(1)+lift(1),d4z(1,outcome_frames(1)+lift(1)),'L')
                subplot(7,1,5)
                text(outcome_frames(1)+lift(1),wz(1,outcome_frames(1)+lift(1)),'L')
                text(outcome_frames(1),max(wz)*0.9,['Lift = ' num2str(outcome_frames(1)+ lift(1))])
            end
        end
        reach = find(d4x(outcome_frames)>-5 & smoothdata(accelvec(t,outcome_frames),'movmean',smoothwin)>accel_thresh/2);
        if isempty(reach)
            reach=nan;
        else
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+reach(1),d4x(1,outcome_frames(1)+reach(1)),'R')
                subplot(7,1,3)
                text(outcome_frames(1)+reach(1),d4z(1,outcome_frames(1)+reach(1)),'R')
                subplot(7,1,5)
                text(outcome_frames(1)+reach(1),wz(1,outcome_frames(1)+reach(1)),'R')
            end
        end
        
        % grasp
        dist3d = d_digdist3d(t,:);
        dmax = find(dist3d==max(dist3d(outcome_frames(1):outcome_frames(1)+500)));
        dmin = find(dist3d==min(dist3d(dmax:outcome_frames(end))));
        grasp_f_all_ind = find(dist3d(dmax:dmin)*100<-1);
        grasp_f = grasp_f_all_ind + dmax;
        if isempty(grasp_f)
            grasp=nan;
        else
            grasp = grasp_f(1)-2;
            if plotfig
                subplot(7,1,3)
                text(grasp(1),d4z(1,grasp(1)),'G')
                subplot(7,1,6)
                text(grasp(1),pz(1,grasp(1)),'G')
            end
        end
        pelletgrasp = find(digpelletdist3d(t,outcome_frames)<=3);
        if isempty(pelletgrasp)
            pelletgrasp=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+pelletgrasp(1),pz(1,outcome_frames(1)+pelletgrasp(1)),'P')
                subplot(7,1,4)
                text(outcome_frames(1)+pelletgrasp(1),pz(1,outcome_frames(1)+pelletgrasp(1)),'P')
            end
        end
        
        digmouth = find((d4x(outcome_frames)>0 & d4x(outcome_frames)<5) & d4y(outcome_frames)>125 & d4z(outcome_frames)>13);
        if isempty(digmouth)
            digmouth=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+digmouth(1),d4z(1,outcome_frames(1)+digmouth(1)),'DM')
            end
        end
        pelletmouth = find((px(outcome_frames)>0 & px(outcome_frames)<4) & py(outcome_frames)>127 & pz(outcome_frames)>15);
        if isempty(pelletmouth)
            pelletmouth=nan;
        else
            if plotfig
                subplot(7,1,3)
                text(outcome_frames(1)+pelletmouth(1),d4z(1,outcome_frames(1)+pelletmouth(1)),'PM')
                subplot(7,1,4)
                text(outcome_frames(1)+pelletmouth(1),pz(1,outcome_frames(1)+pelletmouth(1)),'PM')
            end
        end
        
        % hand open
        % interdigit distance
        dist3d = d_digdist3d(t,:);
        hopen_f = find(dist3d(outcome_frames)*100>4);
        if ~isempty(hopen_f)
            hopen = hopen_f(1)+outcome_frames;
            if plotfig
                subplot(7,1,1)
                text(outcome_frames(1)+hopen(1),d4x(1,outcome_frames(1)+hopen(1)),'H')
            end
        else
            hopen = nan;
        end
        
        
        %supination
        % hand angle
        slope = (d2y - d4y) ./ (d2z - d4z);
        angle = atand(slope);
        angle(angle>0)=90-angle(angle>0);
        angle(angle<0)=90-angle(angle<0);
        ang = smoothdata(angle,'movmean',smoothwin);
        d_ang=[0];
        for a=1:length(ang)
            if a>1
                d_ang(a) = ang(a) - ang(a-1);
            else
                d_ang(a)=0;
            end
        end
        if ~isnan(grasp)
        dmin = find(d_ang==min(d_ang(grasp(1)-10:grasp(1)+50)));
        end
        if ~isempty(dmin)
            supinate = dmin;
            if plotfig
                subplot(7,1,5)
                text(supinate(1),wz(1,supinate),'S')
            end
        else
            supinate = nan;
        end
                
        if plotfig
            figure; set(gcf,'position',[1800 110 1500 500]);
            hold on;
            plot(outcome_frames,dist3d(outcome_frames).*100);
            title('Change in Inter-digit Distance,grasp=star')
            xlabel('frame')
            plot(grasp(1),dist3d(grasp(1))*100,'m*')
            text(outcome_frames(1)+50,max(dist3d(outcome_frames))*100*0.9,['grasp = ' num2str(grasp(1))]);
        end
    end
    
    % if grasping pellet for 200ms, dig at mouth for 50ms, and pellet at
    % mouth for 50ms, then success, otherwise fail
    if length(pelletgrasp)>=50 && length(digmouth)>=25 && length(pelletmouth)>=25
        trial_outcome(t) = 1
    else
        trial_outcome(t) = 0
    end
    
%     if t==110
%         pause()
%     end
    
    % create lift trial index and a flag if its a lift trial (lift_flag)
    % which finds trial_start with lifts occuring prior to length of
    % outcome_frames range
    if ~isempty(lift) & ~isnan(lift)
        lift_t(t,1:length(lift)) = lift;
        if lift_t(t,1)<outcome_frames(end)
            lift_flag(t,1)=1;
        end
    else
        lift_t(t,1)=nan;
        lift_flag(t,1)=0;
    end
    
    % reach
    if ~isempty(reach) & ~isnan(reach)
        reach_t(t,1:length(reach)) = reach;
        if reach_t(t,1)<outcome_frames(end)
            reach_flag(t,1)=1;
        end
    else
        reach_t(t,1)=nan;
        reach_flag(t,1)=0;
    end
    
    % hand open
    if ~isempty(hopen) & ~isnan(hopen)
        hopen_t(t,1:length(hopen)) = hopen;
        if hopen_t(t,1)<outcome_frames(end)
            hopen_flag(t,1)=1;
        end
    else
        hopen_t(t,1)=nan;
        hopen_flag(t,1)=0;
    end
    
    % grasp
    if ~isempty(grasp) & ~isnan(grasp)
        grasp_t(t,1:length(grasp)) = grasp;
        if grasp_t(t,1)<outcome_frames(end)
            grasp_flag(t,1)=1;
        end
    else
        grasp_t(t,1)=nan;
        grasp_flag(t,1)=0;
    end
    
    % pelletgrasp
    if ~isempty(pelletgrasp) & ~isnan(pelletgrasp)
        pelletgrasp_t(t,1:length(pelletgrasp)) = pelletgrasp;
        if pelletgrasp_t(t,1)<outcome_frames(end)
            pelletgrasp_flag(t,1)=1;
        end
    else
        pelletgrasp_t(t,1)=nan;
        pelletgrasp_flag(t,1)=0;
    end
    
    % supinate
    if ~isempty(supinate) & ~isnan(supinate)
        supinate_t(t,1:length(supinate)) = supinate;
        if supinate_t(t,1)<outcome_frames(end)
            supinate_flag(t,1)=1;
        end
    else
        supinate_t(t,1)=nan;
        supinate_flag(t,1)=0;
    end
    
    % pelletmouth
    if ~isempty(pelletmouth) & ~isnan(pelletmouth)
        pelletmouth_t(t,1:length(pelletmouth)) = pelletmouth;
        if pelletmouth_t(t,1)<outcome_frames(end)
            pelletmouth_flag(t,1)=1;
        end
    else
        pelletmouth_t(t,1)=nan;
        pelletmouth_flag(t,1)=0;
    end
    
    % digmouth
    if ~isempty(digmouth) & ~isnan(digmouth)
        digmouth_t(t,1:length(digmouth)) = digmouth;
        if digmouth_t(t,1)<outcome_frames(end)
            digmouth_flag(t,1)=1;
        end
    else
        digmouth_t(t,1)=nan;
        digmouth_flag(t,1)=0;
    end
    
    % Reach velocity
    if ~isempty(reach) & ~isnan(reach)
        if ~isempty(pelletgrasp) & ~isnan(pelletgrasp)
            % distance between pellet and digit / time frome lift to grasp
            digpelletdist3d_curr = norm([d4x(outcome_frames(lift(1))) d4y(outcome_frames(lift(1))) d4z(outcome_frames(lift(1)))] - [d4x(outcome_frames(pelletgrasp(1))) d4y(outcome_frames(pelletgrasp(1))) d4z(outcome_frames(pelletgrasp(1)))]);
            reach_vel(t,1) = digpelletdist3d_curr / (outcome_frames(pelletgrasp(1)) - outcome_frames(lift(1)));
        else
            % distance between digit and digit / time frome lift to reach
            digpelletdist3d_curr = norm([d4x(outcome_frames(lift(1))) d4y(outcome_frames(lift(1))) d4z(outcome_frames(lift(1)))] - [d4x(outcome_frames(reach(1))) d4y(outcome_frames(reach(1))) d4z(outcome_frames(reach(1)))]);
            reach_vel(t,1) = digpelletdist3d_curr / (outcome_frames(reach(1)) - outcome_frames(lift(1)));
        end
    else
        reach_vel(t,1)=nan;
    end
    
    % reaction time
    if ~isempty(lift) & ~isnan(lift)
        reaction_time(t,1) = outcome_frames(lift(1)) - outcome_frames(1);
    else
        reaction_time(t,1)=nan;
    end
    
    % trial outcome array contains some useful variables
    trial_outcome_array(t,:) = [trial_outcome(t) reach_flag(t) lift_flag(t)...
        reach_vel(t) reaction_time(t) outcome_frames(1)+lift_t(1) ...
        outcome_frames(1)+reach_t(1) outcome_frames(1)+grasp(1) ...
        outcome_frames(1)+pelletgrasp(1) outcome_frames(1)+digmouth(1) ...
        outcome_frames(1)+pelletmouth(1)];
  
    if plotfig
        % show grasping dynamics, inter digit distance and change
%         pause();
        close all
    end
end
toc;

tic;
% second loop - compare trial traj to mean of all success control trials
% mean of all success control trial
contrial = con_trial(find(con_trial<=trial_range(end))); % use trial_range
con_succ_trial = contrial(find(trial_outcome(contrial)==1));
d1all = mean(smoothdata(squeeze(traj_2(con_succ_trial,1,:))','movmean',lift_smoothwin),2);
d2all = mean(smoothdata(squeeze(traj_2(con_succ_trial,2,:))','movmean',lift_smoothwin),2);
d3all = mean(smoothdata(squeeze(traj_2(con_succ_trial,3,:))','movmean',lift_smoothwin),2);
for t=1:length(trial_range)
    if strcmp(neural_dir,[rootdir '\jcr82\ephys\jcr82_20210101_4400um_g0']) && t>102
        vid_offset=-3;
    else
        vid_offset=0;
    end
    
    % dig2
    d2x = smoothdata(squeeze(traj_1(t+vid_offset,1,:))','movmean',lift_smoothwin);
    d2y = smoothdata(squeeze(traj_1(t+vid_offset,2,:))','movmean',lift_smoothwin);
    d2z = smoothdata(squeeze(traj_1(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % dig4
    d4x = smoothdata(squeeze(traj_2(t+vid_offset,1,:))','movmean',lift_smoothwin);
    d4y = smoothdata(squeeze(traj_2(t+vid_offset,2,:))','movmean',lift_smoothwin);
    d4z = smoothdata(squeeze(traj_2(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % wrist
    wx = smoothdata(squeeze(traj_3(t+vid_offset,1,:))','movmean',lift_smoothwin);
    wy = smoothdata(squeeze(traj_3(t+vid_offset,2,:))','movmean',lift_smoothwin);
    wz = smoothdata(squeeze(traj_3(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % pellet
    px = smoothdata(squeeze(traj_4(t+vid_offset,1,:))','movmean',lift_smoothwin);
    py = smoothdata(squeeze(traj_4(t+vid_offset,2,:))','movmean',lift_smoothwin);
    pz = smoothdata(squeeze(traj_4(t+vid_offset,3,:))','movmean',lift_smoothwin);
    
    % trajectory difference, linear distance in each dimension and 3d
    % current trial from each dimension
    d1 = d4x; % x
    d2 = d4y;
    d3 = d4z;
    % distance from control
    dist1d_1all=[]; dist1d_2all=[]; dist1d_3all=[]; dist3d_all=[];
    for d=1:length(frames)
        dist1d_1all(d) =  norm([d1(d)] - [d1all(d)]);
    end
    for d=1:length(frames)
        dist1d_2all(d) = norm([d2(d)] - [d2all(d)]);
    end
    for d=1:length(frames)
        dist1d_3all(d) = norm([d3(d)] - [d3all(d)]);
    end
    for d=1:length(frames)
        dist3d_all(d) = norm([d1(d) d2(d) d3(d)] - [d1all(d) d2all(d) d3all(d)]);
    end
    if plotfig
        subplot(7,1,6); plot(dist3d_all(1:length(frames)),'k-'); hold on;
        ylabel('dist 3d')
    end
    
    dist1d_x(t,1:length(frames))=dist1d_1all(1:length(frames));
    dist1d_y(t,1:length(frames))=dist1d_2all(1:length(frames));
    dist1d_z(t,1:length(frames))=dist1d_3all(1:length(frames));
    dist3d(t,1:length(frames))=dist3d_all(1:length(frames));
    
    d2x_i(t,:) = d2x; d2y_i(t,:) = d2y; d2z_i(t,:) = d2z;
    d4x_i(t,:) = d4x; d4y_i(t,:) = d4y; d4z_i(t,:) = d4z;
    wx_i(t,:) = wx; wy_i(t,:) = wy; wz_i(t,:) = wz;
    px_i(t,:) = px; py_i(t,:) = py; pz_i(t,:) = pz;
    
end
toc;

% create lists of trial types
LT_trial= sort([LT_40hz ; LT_10hz ; LT_4hz ; LT_single]);% list of all light trials
LT_lift_trial = lift_flag(sort(LT_trial));% light trials with lift
LT40_lift_trial = lift_flag(sort(LT_40hz));% 40hz light trials with lift
LT10_lift_trial = lift_flag(sort(LT_10hz));% 10hz light trials with lift
LT4_lift_trial = lift_flag(sort(LT_4hz));% 4hz light trials with lift
LTsingle_lift_trial = lift_flag(sort(LT_single));% 4hz light trials with lift
con_lift_trial = lift_flag(contrial);% control trials with lift

LT_reach_trial = reach_flag(sort(LT_trial));% light trials with reach
LT40_reach_trial = reach_flag(sort(LT_40hz));% 40hz light trial with reach
LT10_reach_trial = reach_flag(sort(LT_10hz));% 10hz light trial with reach
LT4_reach_trial = reach_flag(sort(LT_4hz));% 4hz light trial with reach
LTsingle_reach_trial = reach_flag(sort(LT_single));% 4hz light trial with reach
con_reach_trial = reach_flag(contrial);% control trial with reach

% % trajectory difference
% traj_diff_x_con = dist1d_x(con_trial(con_reach_trial==1),:);
% traj_diff_y_con = dist1d_y(con_trial(con_reach_trial==1),:);
% traj_diff_z_con = dist1d_z(con_trial(con_reach_trial==1),:);
% traj_diff_3d_con = dist3d(con_trial(con_reach_trial==1),:);
%
% traj_diff_x_LT40 = dist1d_x(LT_40hz(LT40_reach_trial==1),:);
% traj_diff_y_LT40 = dist1d_y(LT_40hz(LT40_reach_trial==1),:);
% traj_diff_z_LT40 = dist1d_z(LT_40hz(LT40_reach_trial==1),:);
% traj_diff_3d_LT40 = dist3d(LT_40hz(LT40_reach_trial==1),:);

% trajectory difference
% up-sample trajectories and trajectory differences from 2ms to 1ms bins
fs=0.4998;
traj_diff_x=[];traj_diff_y=[];traj_diff_z=[]; traj_diff_3d=[];
traj_d2x=[]; traj_d2y=[]; traj_d2z=[]; traj_d4x=[]; traj_d4y=[]; traj_d4z=[];
traj_wx=[]; traj_wy=[]; traj_wz=[]; traj_px=[]; traj_py=[]; traj_pz=[];
traj_digdist3d=[]; traj_d_digdist3d=[]; traj_digvel=[]; traj_digaccel=[];
traj_digpelletdist3d=[]; traj_wvel=[]; traj_waccel=[]; 
traj_accelvec=[];
    
for i=1:size(d2x_i,1)
    
    traj_d2x(i,:) = interp1(1:size(d2x_i,2),d2x_i(i,:),1:fs:2000);
    traj_d2y(i,:) = interp1(1:size(d2y_i,2),d2y_i(i,:),1:fs:2000);
    traj_d2z(i,:) = interp1(1:size(d2z_i,2),d2z_i(i,:),1:fs:2000);
    
    traj_d4x(i,:) = interp1(1:size(d4x_i,2),d2x_i(i,:),1:fs:2000);
    traj_d4y(i,:) = interp1(1:size(d4y_i,2),d4y_i(i,:),1:fs:2000);
    traj_d4z(i,:) = interp1(1:size(d4z_i,2),d4z_i(i,:),1:fs:2000);
    
    traj_wx(i,:) = interp1(1:size(wx_i,2),wx_i(i,:),1:fs:2000);
    traj_wy(i,:) = interp1(1:size(wy_i,2),wy_i(i,:),1:fs:2000);
    traj_wz(i,:) = interp1(1:size(wz_i,2),wz_i(i,:),1:fs:2000);
    
    traj_px(i,:) = interp1(1:size(px_i,2),px_i(i,:),1:fs:2000);
    traj_py(i,:) = interp1(1:size(py_i,2),py_i(i,:),1:fs:2000);
    traj_pz(i,:) = interp1(1:size(pz_i,2),pz_i(i,:),1:fs:2000);
    
    traj_diff_x(i,:) = interp1(1:size(dist1d_x,2),dist1d_x(i,:),1:fs:2000);
    traj_diff_y(i,:) = interp1(1:size(dist1d_y,2),dist1d_y(i,:),1:fs:2000);
    traj_diff_z(i,:) = interp1(1:size(dist1d_z,2),dist1d_z(i,:),1:fs:2000);
    traj_diff_3d(i,:) = interp1(1:size(dist3d,2),dist3d(i,:),1:fs:2000);
    
    traj_digdist3d(i,:) = interp1(1:size(digdist3d,2),digdist3d(i,:),1:fs:2000);
    traj_d_digdist3d(i,:) = interp1(1:size(d_digdist3d,2),d_digdist3d(i,:),1:fs:2000);
    
    traj_digvel(i,:) = interp1(1:size(digvel,2),digvel(i,:),1:fs:2000);
    traj_digaccel(i,:) = interp1(1:size(digaccel,2),digaccel(i,:),1:fs:2000);
    
    traj_digpelletdist3d(i,:) = interp1(1:size(digpelletdist3d,2),digpelletdist3d(i,:),1:fs:2000);
    
    traj_wvel(i,:) = interp1(1:size(wvel,2),wvel(i,:),1:fs:2000);
    traj_waccel(i,:) = interp1(1:size(waccel,2),waccel(i,:),1:fs:2000);
    
    traj_accelvec(i,:) = interp1(1:size(accelvec,2),accelvec(i,:),1:fs:2000);
    
end

% success rate, within outcome frames window - for reach trial_start
success_rate_con = sum(trial_outcome(contrial))/length(contrial);
success_rate_LT40 = sum(trial_outcome(LT_40hz))/length(LT_40hz);
success_rate_LT10 = sum(trial_outcome(LT_10hz))/length(LT_10hz);

% reach rate, within outcome frames window
reach_rate_con = sum(con_reach_trial)/length(contrial);
reach_rate_LT40 = sum(LT40_reach_trial)/length(LT_40hz);
reach_rate_LT10 = sum(LT10_reach_trial)/length(LT_10hz);

% reach velocity
reach_vel_con = nanmean(reach_vel(contrial(con_reach_trial==1)));
reach_vel_LT40 = nanmean(reach_vel(LT_40hz(LT40_reach_trial==1)));
reach_vel_LT10 = nanmean(reach_vel(LT_10hz(LT10_reach_trial==1)));

% reaction time
reaction_time_con = nanmean(reaction_time(contrial(con_lift_trial==1)));
reaction_time_LT40 = nanmean(reaction_time(LT_40hz(LT40_lift_trial==1)));
reaction_time_LT10 = nanmean(reaction_time(LT_10hz(LT10_lift_trial==1)));


% trial_outcome_struct
% Multiply all timestamp values (in frame number) by 2 to convert to
% timestamp in 1ms bins values. eg. 75 bin = 150 ms
trial_outcome_struct = struct('trial_outcome',trial_outcome',...
    'reach_flag',reach_flag,...
    'lift_flag',lift_flag,...
    'reach_vel',reach_vel,...
    'reaction_time',reaction_time.*2,...
    'lift',lift_t.*2,...
    'reach',reach_t.*2,...
    'grasp',grasp_t.*2,...
    'pelletgrasp',pelletgrasp_t.*2,...
    'digmouth',digmouth_t.*2,...
    'pelletmouth',pelletmouth_t.*2,...
    'handopen',hopen_t*2,...
    'supinate',supinate_t*2);

trajectory_struct = struct('dig2_x',traj_d2x,'dig2_y',traj_d2y,'dig2_z',traj_d2z,...
    'dig4_x',traj_d4x,'dig4_y',traj_d4y,'dig4_z',traj_d4z,...
    'wrist_x',traj_wx,'wrist_y',traj_wy,'wrist_z',traj_wz,...
    'pellet_x',traj_px,'pellet_y',traj_py,'pellet_z',traj_pz,...
    'dist_from_con_x',traj_diff_x,'dist_from_con_y',traj_diff_y,...
    'dist_from_con_z',traj_diff_z,'dist_from_con_3d',traj_diff_3d,...
    'dig2_4_dist',traj_digdist3d,'dig2_4_dist_diff', traj_d_digdist3d,...
    'dig4_vel', traj_digvel, 'dig4_accel', traj_digaccel,'dig4_pellet_dist', traj_digpelletdist3d,...
    'wrist_vel', traj_wvel, 'wrist_accel', traj_waccel,'accel_vector', traj_accelvec);

tic;
% save variable arrays and mean values
if saveData
    % note some variables not listed here may be previously saved in
    % 'create_trial_structure_jcr.m'
    save([save_dir '\event_ind.mat'], ...
        'con_lift_trial', 'con_reach_trial',...
        'LT_trial','LT_lift_trial','LT_reach_trial',...
        'LT40_lift_trial','LT10_reach_trial',...
        'LT4_lift_trial','LT4_reach_trial',...
        'LTsingle_lift_trial','LTsingle_reach_trial',...
        'reach_vel_con','reach_vel_LT40', ...
        'reach_rate_con','reach_rate_LT40', ...
        'reaction_time_con','reaction_time_LT40', ...
        'success_rate_con','success_rate_LT40',...
        'trajectory_struct','trial_outcome_struct','-append');
end
toc;

disp([neural_dir ' - create_trial_outcome complete'])
