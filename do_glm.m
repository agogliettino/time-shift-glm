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
param_all = cell(numCell,1);
variables_all = cell(numCell,1);
AllVarShifts = cell(numCell,1);
AllPShifts = cell(numCell,1);
AllParamVar = cell(numCell,1);
AllVar_TotalLLH_1 = cell(numCell,1);
AllForwardPval = cell(numCell,1);
AllForwardFinal_Pval = cell(numCell,1);
AllvarAllmodelFits = cell(numCell,1);
AmatAll = cell(numCell,1);
Allposgrid = cell(numCell,1);
Allhdgrid = cell(numCell,1);
Allspeedgrid = cell(numCell,1);
Allthetagrid = cell(numCell,1);
AllposVec = cell(numCell,1);
AllhdVec = cell(numCell,1);
AllspdVec = cell(numCell,1);
AllthetaVec = cell(numCell,1);
tuning_curves_all = cell(numCell,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the for loop
encodingCells = encodingCells;
%for n = encodingCells'  
for n = 4
    tic
    %%%%%%%%%%%%%%%%%%%%%%%% LOAD ALL FILES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find directory
    session_dir = Session{n};
    
    session_dir_split =  strsplit(session_dir,'\');
    session_dir_short = session_dir_split{end};
    animalName = session_dir_split{end-2};
    
    % NOTE: this will have to be changed
    %prefix = 'G:\My Drive\All_Cells\'; %this is for the lab machine (windows)
    prefix = '/Volumes/GoogleDrive/My Drive/All_Cells/'; %from my machine (mac)
    
    %find pos file
    posfile = strcat(prefix,animalName,'_',session_dir_short,'_pos.mat');
    
    %find spike file
    spikefile = strcat(prefix,animalName,'_',session_dir_short,'_T',num2str(Tetrode{n}),'C',num2str(Unit{n}),'.mat');
    
    % find eeg file
    eegfile = strcat(prefix,animalName,'_',session_dir_short,'_eeg.mat');
    load(eegfile);

    %%%%%%%%%%%%%%%%%%%%%%%% COLLECT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [param, tuning_curves, variables,VarShifts,P_ShiftVar,paramVar,VarTotal_LLH_1,ForwardPval,ForwardFinal_Pval,varAllModelFits,A,posgrid,hdgrid,speedgrid,thetagrid,posVec,hdVec,spdVec,thetaVec] = create_glm(posfile,spikefile,filt_eeg,sampleRate,BoxSize{n});
    param_all{n} = param;
    tuning_curves_all{n} = tuning_curves;
    variables_all{n} = variables;
    AllVarShifts{n} = VarShifts;
    AllPShifts{n} = P_ShiftVar;
    AllParamVar{n} = paramVar;
    AllVar_TotalLLH_1{n} = VarTotal_LLH_1;
    AllForwardPval{n} = ForwardPval;
    AllForwardFinal_Pval{n} = ForwardFinal_Pval;
    AllvarAllmodelFits{n} = varAllModelFits;
    
    keyboard
    
    toc
    %save('do_glm_output_112817.mat','param_all','variables_all','AllVarShifts','AllPShifts','AllParamVar','AllVar_TotalLLH_1','AllForwardFinal_Pval','AllvarAllmodelFits','tuning_curves_all');
     
        
    
end




return