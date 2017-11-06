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
LLH_all = cell(numCell,1);
param_all = cell(numCell,1);
variables_all = cell(numCell,1);
shift_values_all = nan(numCell,3);
shift_pval_all = nan(numCell,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the for loop
encodingCells = encodingCells(encodingCells >= 542);
for n = encodingCells'
    tic
    %%%%%%%%%%%%%%%%%%%%%%%% LOAD ALL FILES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find directory
    session_dir = Session{n};
    
    session_dir_split =  strsplit(session_dir,'\');
    session_dir_short = session_dir_split{end};
    animalName = session_dir_split{end-2};
    
    % NOTE: this will have to be changed
    prefix = '/Users/kiah/Desktop/LNP_Data/';
    
    %find pos file
    posfile = strcat(prefix,animalName,'_',session_dir_short,'_pos.mat');
    
    %find spike file
    spikefile = strcat(prefix,animalName,'_',session_dir_short,'_T',num2str(Tetrode{n}),'C',num2str(Unit{n}),'.mat');
    
    % find eeg file
    eegfile = strcat(prefix,animalName,'_',session_dir_short,'_eeg.mat');
    load(eegfile);

    %%%%%%%%%%%%%%%%%%%%%%%% COLLECT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [LLH, param, variables, shift_values, shift_pval] = create_glm(posfile,spikefile,filt_eeg,sampleRate,BoxSize{n});
    LLH_all{n} = LLH;
    param_all{n} = param;
    variables_all{n} = variables;
    shift_values_all(n,:) = shift_values;
    shift_pval_all(n,:) = shift_pval;
    toc
    save('time_shift_output_542on.mat','LLH_all','param_all','variables_all','shift_values_all','shift_pval_all')
    
end




return