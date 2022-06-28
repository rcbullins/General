% jcr_MLstart.m
%%
% Run this at first before doing analysis in *_analysis_log*.m

% Global paths

%%%%%%
% select one of the two options for GlobalDataDir
% GlobalDataDir = 'Z:\Data\Jeremy\Intracellular';
GlobalDataDir = 'D:\localData';
GlobalResultsDir = 'D:\localResults';
GlobalProgramsDir = 'D:\matlab';

% Add paths to Programs
addpath(genpath(GlobalProgramsDir))  % recursively add subdirectories
% Set current directory to Results directory
cd(GlobalDataDir)

%
set(0,'DefaultFigureRenderer','painters');
set(0,'DefaultTextInterpreter','none');
format longG

