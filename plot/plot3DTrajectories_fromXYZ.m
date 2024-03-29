function [] = plot3DTrajectories_fromXYZ(x,y,z,pellet,liftFrames, handOpenFrames,grabFrames,atMouthFrames,varargin)
% PURPOSE
%   Plot all numper of sample trials in grey, and mean in red overlayed.
%   Plot mean jaaba behavior for all trials along mean trajectory - this is
%   chosen based off average frame away from lift. Trajectory starts at
%   lift.
% 
%   Will either plot sample trials with mean overlayed, or only mean by
%   itself.
% INPUT
%   XYZ information of object
%   pellet x,y,z information
%   PlotFrameAfter   how many frames to plot after lift (default 1000)
%   color      (default = red)
%   liftFrames  which frame lift occured in for each trial (1 x num trial)
%               beginning of plotting
%   graspFrames
%   handOpenFrames
%   atMouthFrames (end of plotting)
%   PlotMeanOnly (0 or 1) only want to plot mean without box, then make it
%   1 (default is 0)
% OUTPUT
%   subplot containing sample trials from the dataset of 3D trajectory.
% HISTORY
%   Reagan Bullins: 8.22.2022
%% Input Parsers
p = inputParser;
addParameter(p,'PlotFrameBefore',400, @isnumeric);
addParameter(p,'PlotFrameAfter',1000, @isnumeric);
addParameter(p,'Color', [1, 0, 0], @isvector)
addParameter(p,'SampleSize',size(x,1), @isnumeric);
addParameter(p,'PlotMeanOnly',0, @isnumeric);
parse(p,varargin{:});
PlotFrameBefore = p.Results.PlotFrameBefore;
PlotFrameAfter = p.Results.PlotFrameAfter;
PlotMeanOnly = p.Results.PlotMeanOnly;
Color = p.Results.Color;
SampleSize = p.Results.SampleSize;
if PlotMeanOnly == 0
    %% plot 3d trajectories
    % Plot sample trajectories: Lift to at mouth
    for i=1:SampleSize
        if isnan(liftFrames(1,i)) || isnan(handOpenFrames(1,i)) || isnan(grabFrames(1,i)) || isnan(atMouthFrames(1,i))
            continue;
        else
            plot3(x(i,liftFrames(1,i):atMouthFrames(1,i)), ...
                y(i,liftFrames(1,i):atMouthFrames(1,i)),...
                z(i,liftFrames(1,i):atMouthFrames(1,i)),'Color', [105/255 105/255 105/255 .1]);
        end
        hold on;
    end
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    axis square;
    box off;
    set(gca,'xticklabel',{[]});
    set(gca,'yticklabel',{[]});
    set(gca,'zticklabel',{[]});
    % Find where lift does not happen and remove from mean
    non_nan_idx_lift = find(~isnan(liftFrames));
    non_nan_idx_hand = find(~isnan(handOpenFrames));
    non_nan_idx_grasp = find(~isnan(grabFrames));
    non_nan_idx_mouth = find(~isnan(atMouthFrames));

    non_nan_idx  = intersect(intersect(intersect(non_nan_idx_lift,non_nan_idx_hand),non_nan_idx_grasp),non_nan_idx_mouth);

    x_tmp_mean = [];
    y_tmp_mean = [];
    z_tmp_mean = [];
    for i = 1:size(x(non_nan_idx),2)
        x_tmp_mean = [x_tmp_mean; x(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
        y_tmp_mean = [y_tmp_mean; y(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
        z_tmp_mean = [z_tmp_mean; z(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
    end
    
else

    % Find where lift does not happen and remove from mean
    non_nan_idx_lift = find(~isnan(liftFrames));
    non_nan_idx_hand = find(~isnan(handOpenFrames));
    non_nan_idx_grasp = find(~isnan(grabFrames));
    non_nan_idx_mouth = find(~isnan(atMouthFrames));
    % Find all non nan idx for all behaviors
    non_nan_idx  = intersect(intersect(intersect(non_nan_idx_lift,non_nan_idx_hand),non_nan_idx_grasp),non_nan_idx_mouth);

    x_tmp_mean = [];
    y_tmp_mean = [];
    z_tmp_mean = [];
    for i = 1:size(x(non_nan_idx),2)
        x_tmp_mean = [x_tmp_mean; x(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
        y_tmp_mean = [y_tmp_mean; y(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
        z_tmp_mean = [z_tmp_mean; z(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):liftFrames(1,non_nan_idx(i))+PlotFrameAfter)];
    end

end

x_mean = mean(x_tmp_mean);
y_mean = mean(y_tmp_mean);
z_mean = mean(z_tmp_mean);

plot3(x_mean,y_mean,z_mean,'r')
hold on;
% Plot dots on where different parts of reach
% Lift
scatter3(x_mean(1),y_mean(1),z_mean(1),'o','filled','MarkerEdgeColor','flat','MarkerFaceColor', Color);
txt_lift = 'Lift';
text(x_mean(1),y_mean(1),z_mean(1),txt_lift);
lift_mean_frame = round(mean(liftFrames(1,non_nan_idx)));

% Hand Open mean frame
handOpen_mean_frame = round(mean(handOpenFrames(1,non_nan_idx)));
handOpen_inRespectToLift = handOpen_mean_frame-lift_mean_frame;

scatter3(x_mean(handOpen_inRespectToLift),...
    y_mean(handOpen_inRespectToLift),...
    z_mean(handOpen_inRespectToLift), 'o','filled');
txt_ho = 'Hand open';
text(x_mean(handOpen_inRespectToLift),y_mean(handOpen_inRespectToLift),z_mean(handOpen_inRespectToLift),txt_ho);

% % Grab mean frame
grab_mean_frame = round(mean(grabFrames(1,non_nan_idx)));
grab_inRespectToLift = grab_mean_frame-lift_mean_frame;

scatter3(x_mean(grab_inRespectToLift),...
    y_mean(grab_inRespectToLift),...
    z_mean(grab_inRespectToLift), 'o','filled');
txt_g = 'Grab';
text(x_mean(grab_inRespectToLift),y_mean(grab_inRespectToLift),z_mean(grab_inRespectToLift),txt_g);


% At Mouth
mouth_mean_frame = round(mean(atMouthFrames(1,non_nan_idx)));
mouth_inRespectToLift = mouth_mean_frame-lift_mean_frame;

scatter3(x_mean(mouth_inRespectToLift),...
    y_mean(mouth_inRespectToLift),...
    z_mean(mouth_inRespectToLift), 'o','filled');
txt_m = 'At mouth';
text(x_mean(mouth_inRespectToLift),y_mean(mouth_inRespectToLift),z_mean(mouth_inRespectToLift),txt_m);

% Pellet Position
x_tmp_pel = [];
y_tmp_pel = [];
z_tmp_pel = [];
for ipel = 1:length(non_nan_idx)
    x_tmp_pel = [x_tmp_pel; pellet.x(non_nan_idx(ipel),handOpenFrames(1,non_nan_idx(ipel)))];
    y_tmp_pel = [y_tmp_pel; pellet.y(non_nan_idx(ipel),handOpenFrames(1,non_nan_idx(ipel)))];
    z_tmp_pel = [z_tmp_pel; pellet.z(non_nan_idx(ipel),handOpenFrames(1,non_nan_idx(ipel)))];
end
x_pel_mean = mean(x_tmp_pel);
y_pel_mean = mean(y_tmp_pel);
z_pel_mean = mean(z_tmp_pel);
scatter3(x_pel_mean, y_pel_mean, z_pel_mean, 'o','filled');
txt_p = 'Pellet';
text(x_pel_mean, y_pel_mean, z_pel_mean,txt_p);


if PlotMeanOnly
   box off;
   set(gca,'xtick',[]);
   set(gca,'ytick',[]);
   set(gca,'ztick',[]);
   axis off;
end


%%%
% x_tmp_mean =[];
%    y_tmp_mean = [];
%    z_tmp_mean=[];
%    interp_pts = 1000;
%     for i = 1:size(x(non_nan_idx),2)
%         this_x = x(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):atMouthFrames(1,non_nan_idx(i)));
%         this_y = y(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):atMouthFrames(1,non_nan_idx(i)));
%         this_z = z(non_nan_idx(i), liftFrames(1,non_nan_idx(i)):atMouthFrames(1,non_nan_idx(i)));
%
%         xi = cumsum([1, abs(diff(this_x))]);
%         xq = min(xi):(max(xi)-min(xi))/interp_pts:max(xi);
%         x_interp = interp1(xi,this_x,xq);
%
%         yi = cumsum([1, abs(diff(this_y))]);
%         yq = min(yi):(max(yi)-min(yi))/interp_pts:max(yi);
%         y_interp = interp1(yi,this_y,yq);
%
%         zi = cumsum([1, abs(diff(this_z))]);
%         zq = min(zi):(max(zi)-min(zi))/interp_pts:max(zi);
%         z_interp = interp1(zi,this_z,zq);
%
%         x_tmp_mean = [x_tmp_mean; x_interp];
%         y_tmp_mean = [y_tmp_mean; y_interp];
%         z_tmp_mean = [z_tmp_mean; z_interp];



%    % Take the mean positions of lift frames to at mouth frames
%     x_mean = mean(x(non_nan_idx,liftFrames(1,non_nan_idx):atMouthFrames(1,non_nan_idx)));
%     y_mean = mean(y(non_nan_idx,liftFrames(1,non_nan_idx):atMouthFrames(1,non_nan_idx)));
%     z_mean = mean(z(non_nan_idx,liftFrames(1,non_nan_idx):atMouthFrames(1,non_nan_idx)));
%     plot3(x_mean, y_mean,z_mean,'Color',Color);
%     hold on;