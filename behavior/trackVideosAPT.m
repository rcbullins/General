function [] = trackVideosAPT(movie_dir,video_perspective,lObj)
% PURPOSE
%   Run tracking over all videos within movie_dir pathway. Will create a
%   'trk' folder within movie directory and store all .trk files for each
%   movie within.
% INPUTS
%   movie_dir:          directory with all movies (string)
%   video_perspective:  'front' or 'side'         (string)
% OUTPUS
%   saved trk files to movie directory.
%% get movie names
tmp_files=dir(movie_dir);
movies={};
for m=1:length(tmp_files)
    curr_fil=tmp_files(m).name;
    
    if ~isempty(strfind(curr_fil,'.avi'))
        if strcmp(video_perspective,'front')
           if ~isempty(strfind(curr_fil,'front')) 
               tmp_count=size(movies,1);
               movies{tmp_count+1,1}=strcat(movie_dir,'\',curr_fil);
               movies{tmp_count+1,2}=strcat(movie_dir,'\trk\','tracked_',video_perspective,'_',curr_fil(end-6:end-4),'.trk')
           end
        else
       
           if ~isempty(strfind(curr_fil,'side')) 
               tmp_count=size(movies,1);
               movies{tmp_count+1,1}=strcat(movie_dir,'\',curr_fil);
               movies{tmp_count+1,2}=strcat(movie_dir,'\trk\','tracked_',video_perspective,'_',curr_fil(end-6:end-4),'.trk')
           end
        end
    end
        
end

% load APT
nfiles=size(movies,1);
mkdir(strcat(movie_dir,'\trk'))
f0=ones([nfiles,1]);
f1=ones([nfiles,1])*3000;

for i=1:20:nfiles %nfiles
    flag=1;
    counter=1;
    while flag
        if ~lObj.tracker.bgTrkIsRunning
            flag=0
            
        end
        pause(10)
        counter=counter+1
    end
    if size(movies,1)>(i+19)
        lObj.tracker.track(movies(i:i+19,1),'trkfiles',movies(i:i+19,2),'f0',f0(1:20),'f1',f1(1:20));
    else
        lObj.tracker.track(movies(i:end,1),'trkfiles',movies(i:end,2),'f0',f0(i:end),'f1',f1(i:end));
    end
    i
    pause(15)
    txt='here'
end
