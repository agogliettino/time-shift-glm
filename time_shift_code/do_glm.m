function [] = do_glm()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load spreadsheet
[~,~,data] = xlsread('Wildtype_Mice_Database_v5.xlsx');

% nighteyes - NE (an animal) - recordings were from deeper layers of MEC.
% excluded from this analyses
withoutNE = 878;
Session = data(2:withoutNE,11); 
Tetrode = data(2:withoutNE,3);
Unit = data(2:withoutNE,4);
BoxSize = data(2:withoutNE,13);
numCell = numel(Tetrode);

load encodingCells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize matrices
testFit_all = cell(numCell,1);
param_all = cell(numCell,1);
variables_all = cell(numCell,1);
ShiftLLHMax_all = cell(numCell,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the for loop

for n = encodingCells'
    %%%%%%%%%%%%%%%%%%%%%%%% LOAD ALL FILES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find directory
    session_dir = Session{n};
    
    session_dir_split =  strsplit(session_dir,'\');
    session_dir_short = session_dir_split{end};
    animalName = session_dir_split{end-2};
    
    % NOTE: this will have to be changed
    prefix = 'G:\My Drive\All_Cells\';
    
    %find pos file
    posfile = strcat(prefix,animalName,'_',session_dir_short,'_pos.mat');
    
    %find spike file
    spikefile = strcat(prefix,animalName,'_',session_dir_short,'_T',num2str(Tetrode{n}),'C',num2str(Unit{n}),'.mat');
    
    % find eeg file
    eegfile = strcat(prefix,animalName,'_',session_dir_short,'_eeg.mat');
    load(eegfile);

    %%%%%%%%%%%%%%%%%%%%%%%% COLLECT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [testFit, param, variables,AllMaxLLHShifts] = create_glm(posfile,spikefile,filt_eeg,sampleRate,BoxSize{n});
    testFit_all{n} = testFit;
    param_all{n} = param;
    variables_all{n} = variables;
    ShiftLLHMax_all{n} = AllMaxLLHShifts;
    
    
end




return