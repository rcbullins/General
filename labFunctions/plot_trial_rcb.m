function plot_trial_jcr(event_dir,trial,t_win,smoothwin,expt_type)
%% load data
load(event_dir)

outcome_frames = 1000:2500; % time range (1ms bins)

%expt_type=exp_type(i);

%%
% expt_type = exp_type(i);

% successful trials
sc=find(trial_outcome_struct.trial_outcome==1);
usc=find(trial_outcome_struct.trial_outcome==0);

rt = find(trial_outcome_struct.reach_flag==1);
[rt_con,pos] = intersect(rt,con_trial);

% unsuccesful controls
[usc_con,pos] = intersect(usc,con_trial);

% succesful controls
[sc_con,pos] = intersect(sc,con_trial);

laser_all=sort([LT_10hz; LT_40hz]);
% successful Laser
[sc_laser,pos] = intersect(sc,laser_all);

% successful LT_40hz
[sc_laser_all] = intersect(sc,LT_40hz);



%%
t_range=sc_laser;
t_range=laser_all;
t_range = con_trial;


for t=1:length(t_range)
    close all;
        
    trial=t_range(t);
    
    if trial_outcome_struct.trial_outcome(trial)
        disp(['trial = ' num2str(trial) ', success'])
    else
        disp(['trial = ' num2str(trial) ', fail'])
    end
    
    
    f1=figure(1);
    set(gcf,'Position',[200 50 1600 500]);
    d2x = trajectory_struct.dig2_x(trial,:)';
    d2y = trajectory_struct.dig2_y(trial,:)';
    d2z = trajectory_struct.dig2_z(trial,:)';
    subplot(4,3,1)
    hold on;
    plot(d2x); % plot all trials
    ylabel('Digit_2')
    title('X')
    subplot(4,3,2)
    hold on;
    plot(d2y); % plot all trials
    title('Y')
    subplot(4,3,3)
    hold on;
    plot(d2z); % plot all trials
    title('Z')
    
    d4x = trajectory_struct.dig4_x(trial,:)';
    d4y = trajectory_struct.dig4_y(trial,:)';
    d4z = trajectory_struct.dig4_z(trial,:)';
    subplot(4,3,4)
    hold on;
    plot(d4x); % plot all trials
    ylabel('Digit_4')
    subplot(4,3,5)
    hold on;
    plot(d4y); % plot all trials
    subplot(4,3,6)
    hold on;
    plot(d4z); % plot all trials
    
    wx = trajectory_struct.wrist_x(trial,:)';
    wy = trajectory_struct.wrist_y(trial,:)';
    wz = trajectory_struct.wrist_z(trial,:)';
    subplot(4,3,7)
    hold on;
    plot(wx); % plot all trials
    ylabel('Wrist')
    subplot(4,3,8)
    hold on;
    plot(wy); % plot all trials
    subplot(4,3,9)
    hold on;
    plot(wz); % plot all trials
    
    px = trajectory_struct.pellet_x(trial,:)';
    py = trajectory_struct.pellet_y(trial,:)';
    pz = trajectory_struct.pellet_z(trial,:)';
    subplot(4,3,10)
    hold on;
    plot(px); % plot all trials
    ylabel('Pellet')
    xlabel('time (ms)')
    subplot(4,3,11)
    hold on;
    plot(py); % plot all trials
    xlabel('time (ms)')
    subplot(4,3,12)
    hold on;
    plot(pz); % plot all trials
    xlabel('time (ms)')
    
    linkaxesInFigure('x')
    subplot(4,3,12)
    xlim([outcome_frames(1) outcome_frames(end)]);
    zoom on;
    
    
    f2=figure(2);
    set(gcf,'Position',[200 400 1600 500]);
    
    % Supination
    colord=copper(length(t_win));
    subplot(1,4,1)
    if 1
        for f=t_win(1):10:t_win(end)
            plot(d2y(f), d2z(f),'r.','markerSize',12);
            hold all;
            plot(d4y(f), d4z(f),'k.','markerSize',12);
            hold all;
            line([d2y(f) d4y(f)],[d2z(f) d4z(f)],'Color',colord(f-t_win(1)+1,:))
            hold all;
        end
    end
    title('Digit-2 and Digit-4 Position')
    xlabel('y position')
    ylabel('z position')
    
    if expt_type==1
        ylim([-11 5]);
    elseif expt_type==2 || expt_type==3
        ylim([0 20]);
    end
    
    % hand angle
    slope = (d2y - d4y) ./ (d2z - d4z);
    angle = atand(slope);
    angle(angle>0)=90-angle(angle>0);
    angle(angle<0)=90-angle(angle<0);
    ang = smoothdata(angle,'movmean',smoothwin);
    d_ang=[0 0 0 0 0];
    for a=1:length(ang)
        if a>1
            d_ang(a,:) = ang(a,:) - ang(a-1,:);
        else
            d_ang(a,:)=0;
        end
    end
    
    subplot(1,4,2)
    cla;
    plot(t_win,ang(t_win)); hold on;
    plot(t_win,d_ang(t_win)*100); hold on;
    
    title('Hand Rotation Angle and Change in Angle')
    xlabel('time (ms)')
    ylim([-90 180]);
    xlim([outcome_frames(1) outcome_frames(end)]);
    
    % inter-digit distance
    dist3d = trajectory_struct.dig2_4_dist(trial,:);
    d_dist3d = trajectory_struct.dig2_4_dist_diff(trial,:);
    % lift
    lift = trial_outcome_struct.lift(trial(1))+outcome_frames(1);
    disp(['lift_time=' num2str(lift) ', lift_frame=' num2str(lift/2)]);
    
    % hand open
    hopen = trial_outcome_struct.handopen(trial,1);
    disp(['hopen_time=' num2str(hopen) ', hopen_frame=' num2str(hopen/2)]);
    
    % grasp
    grasp = trial_outcome_struct.grasp(trial(1));
    disp(['grasp_time=' num2str(grasp) ', grasp_frame=' num2str(grasp/2)]);
    
    % supinate
    supinate = trial_outcome_struct.supinate(trial(1));
    disp(['supinate_time=' num2str(supinate) ', supinate_frame=' num2str(supinate/2)]);
    
    % pellet at mouth
    pm = trial_outcome_struct.pelletmouth(trial(1))+outcome_frames(1);
    disp(['atmouth_time=' num2str(pm) ', atmouth_frame=' num2str(pm/2)]);
    
    subplot(1,4,4);
    cla
    plot3(d4x,d4y,d4z); hold on;
    title('Digit 4 Position')
    
    subplot(1,4,3)
    cla
    hold on;
    plot(outcome_frames,d_dist3d(1,outcome_frames).*100);
    plot(outcome_frames,dist3d(1,outcome_frames));
    title('Inter-digit Distance and Change in Distance')
    xlabel('time (ms)')
    
    
    try
        
        subplot(1,4,4);
        plot3(d4x(lift), d4y(lift), d4z(lift),'*','color',[1 0 1])
        plot3(d4x(hopen), d4y(hopen), d4z(hopen),'*','color',[1 0 0.8])
        plot3(d4x(grasp), d4y(grasp), d4z(grasp),'*','color',[1 0 0.6])
        plot3(d4x(supinate), d4y(supinate), d4z(supinate),'*','color',[1 0 0.4])
        plot3(d4x(pm), d4y(pm), d4z(pm),'*','color',[1 0 0.2])
        
        subplot(1,4,3)
        plot(lift,d_dist3d(lift)*100,'*','color',[1 0 1])
        plot(hopen,d_dist3d(hopen)*100,'*','color',[1 0 0.8])
        plot(grasp,d_dist3d(grasp)*100,'*','color',[1 0 0.6])
        plot(supinate,d_dist3d(supinate)*100,'*','color',[1 0 0.4])
        plot(pm,d_dist3d(pm)*100,'*','color',[1 0 0.2])
        
        tloc = max(get(gca,'ylim'))*0.9;
        text(lift,tloc,'L','color',[0 0 1])
        text(hopen,tloc,'H','color',[0 0 1])
        text(grasp,tloc,'G','color',[0 0 1])
        text(supinate,tloc,'S','color',[0 0 1])
        text(pm,tloc,'M','color',[0 0 1])
        
        subplot(1,4,2)
        plot(lift,d_ang(lift)*100,'*','color',[1 0 1])
        plot(hopen,d_ang(hopen)*100,'*','color',[1 0 0.8])
        plot(grasp,d_ang(grasp)*100,'*','color',[1 0 0.6])
        plot(supinate,d_ang(supinate)*100,'*','color',[1 0 0.4])
        plot(pm,d_ang(pm)*100,'*','color',[1 0 0.2])
        
        tloc = max(get(gca,'ylim'))*0.9;
        text(lift,tloc,'L','color',[0 0 1])
        text(hopen,tloc,'H','color',[0 0 1])
        text(grasp,tloc,'G','color',[0 0 1])
        text(supinate,tloc,'S','color',[0 0 1])
        text(pm,tloc,'M','color',[0 0 1])
        
        figure(f1);
        subplot(4,3,1)
        plot(lift,d2x(lift),'*','color',[1 0 1])
        plot(hopen,d2x(hopen),'*','color',[1 0 0.8])
        plot(grasp,d2x(grasp),'*','color',[1 0 0.6])
        plot(supinate,d2x(supinate),'*','color',[1 0 0.4])
        plot(pm,d2x(pm),'*','color',[1 0 0.2])
        subplot(4,3,4)
        plot(lift,d4x(lift),'*','color',[1 0 1])
        plot(hopen,d4x(hopen),'*','color',[1 0 0.8])
        plot(grasp,d4x(grasp),'*','color',[1 0 0.6])
        plot(supinate,d4x(supinate),'*','color',[1 0 0.4])
        plot(pm,d4x(pm),'*','color',[1 0 0.2])
        subplot(4,3,7)
        plot(lift,wx(lift),'*','color',[1 0 1])
        plot(hopen,wx(hopen),'*','color',[1 0 0.8])
        plot(grasp,wx(grasp),'*','color',[1 0 0.6])
        plot(supinate,wx(supinate),'*','color',[1 0 0.4])
        plot(pm,wx(pm),'*','color',[1 0 0.2])
        subplot(4,3,10)
        plot(lift,px(lift),'*','color',[1 0 1])
        plot(hopen,px(hopen),'*','color',[1 0 0.8])
        plot(grasp,px(grasp),'*','color',[1 0 0.6])
        plot(supinate,px(supinate),'*','color',[1 0 0.4])
        plot(pm,px(pm),'*','color',[1 0 0.2])
        subplot(4,3,3)
        plot(lift,d2z(lift),'*','color',[1 0 1])
        plot(hopen,d2z(hopen),'*','color',[1 0 0.8])
        plot(grasp,d2z(grasp),'*','color',[1 0 0.6])
        plot(supinate,d2z(supinate),'*','color',[1 0 0.4])
        plot(pm,d2z(pm),'*','color',[1 0 0.2])
        subplot(4,3,6)
        plot(lift,d4z(lift),'*','color',[1 0 1])
        plot(hopen,d4z(hopen),'*','color',[1 0 0.8])
        plot(grasp,d4z(grasp),'*','color',[1 0 0.6])
        plot(supinate,d4z(supinate),'*','color',[1 0 0.4])
        plot(pm,d4z(pm),'*','color',[1 0 0.2])
        subplot(4,3,9)
        plot(lift,wz(lift),'*','color',[1 0 1])
        plot(hopen,wz(hopen),'*','color',[1 0 0.8])
        plot(grasp,wz(grasp),'*','color',[1 0 0.6])
        plot(supinate,wz(supinate),'*','color',[1 0 0.4])
        plot(pm,wz(pm),'*','color',[1 0 0.2])
        subplot(4,3,12)
        plot(lift,pz(lift),'*','color',[1 0 1])
        plot(hopen,pz(hopen),'*','color',[1 0 0.8])
        plot(grasp,pz(grasp),'*','color',[1 0 0.6])
        plot(supinate,pz(supinate),'*','color',[1 0 0.4])
        plot(pm,pz(pm),'*','color',[1 0 0.2])
    catch
    end
    
    pause();
end

%%
figure();
plot(LT_10hz,ones(1,length(LT_10hz)),'*'); hold on; plot(LT_40hz,ones(1,length(LT_40hz)),'*'); set(gca,'xlim',[0 180],'ylim',[0 3])

g1=sum(trial_outcome_struct.reach_flag(con_trial(1:16)))/16*100
g2=sum(trial_outcome_struct.reach_flag(con_trial(17:32)))/16*100
g3=sum(trial_outcome_struct.reach_flag(con_trial(33:48)))/16*100
g4=sum(trial_outcome_struct.reach_flag(con_trial(49:65)))/16*100
g5=sum(trial_outcome_struct.reach_flag(con_trial(66:81)))/16*100
g6=sum(trial_outcome_struct.reach_flag(con_trial(82:97)))/16*100
g7=sum(trial_outcome_struct.reach_flag(con_trial(98:114)))/16*100

h1=sum(trial_outcome_struct.trial_outcome(con_trial(1:16)*2))/16*100
h2=sum(trial_outcome_struct.trial_outcome(con_trial(17:32)))/16*100
h3=sum(trial_outcome_struct.trial_outcome(con_trial(33:48)))/16*100
h4=sum(trial_outcome_struct.trial_outcome(con_trial(49:65)))/16*100
h5=sum(trial_outcome_struct.trial_outcome(con_trial(66:81)))/16*100
h6=sum(trial_outcome_struct.trial_outcome(con_trial(82:97)))/16*100
h7=sum(trial_outcome_struct.trial_outcome(con_trial(98:114)))/16*100

figure();
subplot(2,1,1);bar(1:6,[h1,h2,h3,h4,h5,h6])
ylabel('success rate')
subplot(2,1,2);bar(1:6,[g1,g2,g3,g4,g5,g6])
ylabel('reach rate')

colors=winter(5);
rm1=mean(trajectory_struct.dig4_z(rt_con(1:10),:),1);
rm2=mean(trajectory_struct.dig4_z(rt_con(11:20),:),1);
rm3=mean(trajectory_struct.dig4_z(rt_con(21:30),:),1);
rm4=mean(trajectory_struct.dig4_z(rt_con(31:40),:),1);
rm5=mean(trajectory_struct.dig4_z(rt_con(41:50),:),1);
figure();
plot(rm1,'color',colors(1,:)); hold on; 
plot(rm2,'color',colors(2,:)); 
plot(rm3,'color',colors(3,:)); 
plot(rm4,'color',colors(4,:)); 
plot(rm5,'color',colors(5,:));
