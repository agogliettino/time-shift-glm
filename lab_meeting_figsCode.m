%code to plot the figures I need for giocomo lab meeting!

clear all; close all; clc %clear the workspace

%% load the data

% load the cells of interest
load('encodingCells.mat');

% load the glm output
load('do_glm_output_112817.mat');

% load excel file with other information

xls_data = xlsread('Wildtype_Mice_Database_v5.xlsx');
boxSizeCell = xls_data(:,13);
cellboxSize100 = find(boxSizeCell == 100);

% to get cells only in the 
encodingCells100 = intersect(encodingCells,cellboxSize100);

%% remove theta from the variables_all, AllPShifts, and AllVarShifts

for i = 1:length(variables_all)
    for j = 1:numel(variables_all{i})
        if variables_all{i}(j) == 4 %if its theta
            variables_all{i}(j) = 0; %set the number to 0
            AllPShifts{i}(j) = []; %make this empty
            AllVarShifts{i}(j) = []; %make this empty
        end
    end
end

% remove zeros from the variables to make it empty

for i = 1:length(variables_all)
    if ismember(0,variables_all{i}) %if there is a zero
        variables_all{i}(find(variables_all{i} == 0)) = []; %make any zeros empty
    end
end



%% find the pos, HD, and spd cells

%these are all the cells which encode these variables
posModCellsAll = [];
pos_shift_all = [];
for k = encodingCells100'
    if ismember(1,variables_all{k}) && AllForwardFinal_Pval{k} < 0.05
        posModCellsAll(end+1) = k; 
       
        pval = AllPShifts{k}{1};
        if pval < 0.05
            pos_shift_all = [pos_shift_all; AllVarShifts{k}(1)];
        else
            pos_shift_all = [pos_shift_all; 0];
        end
        
    end
end

%hd
hdModCellsAll = [];
hd_shift_all = [];
for k = encodingCells100'
    if ismember(2,variables_all{k}) && AllForwardFinal_Pval{k} <.05
        hdModCellsAll(end+1) = k;
        if variables_all{k} == 2
            pval = AllPShifts{k}{1};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(1)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 2 && isequal(variables_all{k}, [1 2])
            pval = AllPShifts{k}{2};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(2)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 2 && (isequal(variables_all{k},[2 3]))
            pval = AllPShifts{k}{1};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(1)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 3 && (isequal(variables_all{k},[1 2 3]))
            pval = AllPShifts{k}{2};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(2)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        end
    
    end
    
end

% spd

spdModCellsAll = [];
spd_shift_all = [];

for k = encodingCells100'
    if ismember(3,variables_all{k}) && AllForwardFinal_Pval{k} < .05
        spdModCellsAll(end+1) = k;
        if variables_all{k} == 3
            pval = AllPShifts{k}{1}; 
            if pval < .05
                spd_shift_all = [spd_shift_all; AllVarShifts{k}];
            else
                spd_shift_all = [spd_shift_all; 0];
            end 
        elseif numel(variables_all{k}) == 2 && (isequal(variables_all{k},[1 3]) || isequal(variables_all{k},[2 3]))
            pval = AllPShifts{k}{2}; %if its 1 3 or 2 3
            if pval <.05
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(2)];
            else
                spd_shift_all = [spd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 2 && isequal(variables_all{k},[3 4])
            pval = AllPShifts{k}{1}; %if its 3 4
            if pval < .05
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(1)];
            else
                spd_shift_all = [spd_shift_all; 0];    
            end
        elseif numel(variables_all{k}) == 3 && isequal(variables_all{k},[1 2 3])
            pval = AllPShifts{k}{3}; %if its 1 2 3
            if pval < .05
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(3)];
            else
                spd_shift_all = [spd_shift_all; 0];    
            end
        elseif numel(variables_all{k}) == 3 && isequal(variables_all{k},[2 3 4])
            pval = AllPShifts{k}{2}; %if its 2 3 4
            if pval < .05
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(2)];
            else
                spd_shift_all = [spd_shift_all; 0];    
            end
        elseif numel(variables_all{k}) == 3 && isequal(variables_all{k},[1 3 4])
            pval = AllPShifts{k}{2};
            if pval < .05 %if its 1 3 4
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(2)];
            else
                spd_shift_all = [spd_shift_all; 0];    
            end
        elseif numel(variables_all{k}) == 4
            pval = AllPShifts{k}{3};
            if pval < .05 %if its 1 2 3 4
                spd_shift_all = [spd_shift_all; AllVarShifts{k}(3)];
            else
                spd_shift_all = [spd_shift_all; 0];    
            end
            
            
        end   
    end
end

%% separate cells which are specifically model encoding

% first, group cells based on the different possible models

posMod = cell(numel(encodingCells100),2); %1
hdMod = cell(numel(encodingCells100),2); %2 
spdMod = cell(numel(encodingCells100),2); %3
posHdMod = cell(numel(encodingCells100),2); %1 2
hdSpdMod = cell(numel(encodingCells100),2); %2 3
posSpdMod = cell(numel(encodingCells100),2); %1 3
posHdSpdMod = cell(numel(encodingCells100),2); %1 2 3

%create arrays for the multiple models - hardcoded, but dependent on how
%the design matrix was created
ph = [1 2];
hs = [2 3];
ps = [1 3];
phs = [1 2 3];

for i = encodingCells100'
    if numel(variables_all{i}) == 1 && AllForwardFinal_Pval{i} <.05 %if there is one variable in the model and the p val is sig do the next line
        if isequal(variables_all{i},1) 
            posMod{i,1} = i;
            if cell2mat(AllPShifts{i}) < .05
                posMod{i,2}(1) = AllVarShifts{i};
            else
                posMod{i,2}(1) = 0;
            end
        elseif isequal(variables_all{i},2) 
            hdMod{i,1} = i;
            if cell2mat(AllPShifts{i}) < .05
               hdMod{i,2}(1) = AllVarShifts{i};
            else
                hdMod{i,2}(1) = 0;
            end
        elseif isequal(variables_all{i},3)
            spdMod{i,1} = i;
            if cell2mat(AllPShifts{i}) < .05
               spdMod{i,2}(1) = AllVarShifts{i};
            else
                spdMod{i,2}(1) = 0;
            end
        end
    elseif numel(variables_all{i}) == 2 && AllForwardFinal_Pval{i} <.05%if there are two variables in the model, run these lines
        if isequal(variables_all{i},ph)
            posHdMod{i,1} = i;
            for m = 1:numel(variables_all{i})
                if cell2mat(AllPShifts{i}(m)) < .05
                   posHdMod{i,2}(m) = AllVarShifts{i}(m);
                else
                    posHdMod{i,2}(m) = 0;
                end
            end    
        elseif isequal(variables_all{i},hs)
            hdSpdMod{i,1} = i;
            for m = 1:numel(variables_all{i})
                if cell2mat(AllPShifts{i}(m)) < .05
                    hdSpdMod{i,2}(m) = AllVarShifts{i}(m);
                else
                    hdSpdMod{i,2}(m) = 0;
                end
            end
        elseif isequal(variables_all{i},ps)
            posSpdMod{i,1} = i;
            for m = 1:numel(variables_all{i})
                if cell2mat(AllPShifts{i}(m)) < .05
                    posSpdMod{i,2}(m) = AllVarShifts{i}(m);
                else
                    posSpdMod{i,2}(m) = 0;
                end
            end  
        end
    elseif numel(variables_all{i}) == 3 && AllForwardFinal_Pval{i} <.05%if there are three variables in the model, run these
        if isequal(variables_all{i},phs) 
            posHdSpdMod{i,1} = i;
            for m = 1:numel(variables_all{i})
                if cell2mat(AllPShifts{i}(m)) < .05
                    posHdSpdMod{i,2}(m) = AllVarShifts{i}(m);
                else
                    posHdSpdMod{i,2}(m) = 0;
                end
            end
        end
    end
end

% convert from cell array to matrix

posMod = cell2mat(posMod);
hdMod = cell2mat(hdMod);
spdMod = cell2mat(spdMod);
posHdMod = cell2mat(posHdMod);
hdSpdMod = cell2mat(hdSpdMod);
posSpdMod = cell2mat(posSpdMod);
posHdSpdMod = cell2mat(posHdSpdMod);

%% plot the different models versus number of cells!

% number of cells in each of the models 
numCellsMod = [numel(posMod(:,1)) numel(hdMod(:,1)) numel(spdMod(:,1)) numel(posHdMod(:,1)) numel(posSpdMod(:,1)) numel(hdSpdMod(:,1)) ...
                numel(posHdSpdMod(:,1))];

%labels indicating the model
modelType = categorical({'P','H','S','PH','PS','HS','PHS'});
modelType = reordercats(modelType,{'P','H','S','PH','PS','HS','PHS'});

%plot the data
figure(1);
bar(modelType,numCellsMod);
ylabel('Number of Cells','FontSize',20);
xlabel('Model Type','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% plot the differnt models vers cells, but detailing the model number of shifts

%the number of each shift for each modeltype
numCellsModByShift = [numel(find(posMod(:,2)>0)) numel(find(posMod(:,2)<0)) numel(find(posMod(:,2)==0)); numel(find(hdMod(:,2)>0)) numel(find(hdMod(:,2)<0)) numel(find(hdMod(:,2)==0));...
                    numel(find(spdMod(:,2)>0)) numel(find(spdMod(:,2)<0)) numel(find(spdMod(:,2)==0)); numel(find(posHdMod(:,2)>0)) numel(find(posHdMod(:,2)<0)) numel(find(posHdMod(:,2)==0));...
                    numel(find(posHdMod(:,3)>0)) numel(find(posHdMod(:,3)<0)) numel(find(posHdMod(:,3)==0)); numel(find(posSpdMod(:,2)>0)) numel(find(posSpdMod(:,2)<0)) numel(find(posSpdMod(:,2)==0));...
                    numel(find(posSpdMod(:,3)>0)) numel(find(posSpdMod(:,3)<0)) numel(find(posSpdMod(:,3)==0)); numel(find(hdSpdMod(:,2)>0)) numel(find(hdSpdMod(:,2)<0)) numel(find(hdSpdMod(:,2)==0));...
                    numel(find(hdSpdMod(:,3)>0)) numel(find(hdSpdMod(:,3)<0)) numel(find(hdSpdMod(:,3)==0)); numel(find(posHdSpdMod(:,2)>0)) numel(find(posHdSpdMod(:,2)<0)) numel(find(posHdSpdMod(:,2)==0));...
                    numel(find(posHdSpdMod(:,3)>0)) numel(find(posHdSpdMod(:,3)<0)) numel(find(posHdSpdMod(:,3)==0)); numel(find(posHdSpdMod(:,4)>0)) numel(find(posHdSpdMod(:,4)<0)) numel(find(posHdSpdMod(:,4)==0))];

%labels
modelTypeByShift = categorical({'P','H','S','PH (Pos)','PH (HD','PS (Pos)', 'PS (Spd)','HS (HD)','HS (Spd)','PHS (Pos)','PHS (HD)','PHS (Spd)'});
modelTypeByShift = reordercats(modelTypeByShift,{'P','H','S','PH (Pos)','PH (HD','PS (Pos)', 'PS (Spd)','HS (HD)','HS (Spd)','PHS (Pos)','PHS (HD)','PHS (Spd)'});          
%plot the data
figure(2);
bar(modelTypeByShift,numCellsModByShift,'stacked');
ylabel('Number of Cells (By Shift Type)','FontSize',20);
xlabel('Model Type','FontSize',20);
legend('Positive Shift','Negative Shift','No Shift');
ax = gca;
ax.FontSize = 20;

%% plot the total pos, Hd, spd models and detail the shift types

numCellAllModsByShift = [numel(find(pos_shift_all>0)) numel(find(pos_shift_all<0)) numel(find(pos_shift_all==0)); numel(find(hd_shift_all>0)) numel(find(hd_shift_all<0)) numel(find(hd_shift_all==0));...
                  numel(find(spd_shift_all>0)) numel(find(spd_shift_all<0)) numel(find(spd_shift_all==0));];

allModsType = categorical({'All Position Encoding Cells','All Head Direction Encoding Cells','All Speed Encoding Cells'});
allModsType = reordercats(allModsType,{'All Position Encoding Cells','All Head Direction Encoding Cells','All Speed Encoding Cells'});
%plot these data
figure(3);
bar(allModsType,numCellAllModsByShift,'stacked');
ylabel('Number of Cells (By Shift Type)','FontSize',20);
xlabel('Pooled Encoding Cells','FontSize',20);
legend('Positive Shift','Negative Shift','No Shift');
ax = gca;
ax.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%NEED TO RUN STATS ON THESE DISTRIBUTIONS TO SEE IF
%%%%%%%%%%%%%%%%%%%%%THEY ARE SIGNIGFICANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%doin the stats

% pos
pBinPosShiftAllPosit = myBinomTest(numel(find(pos_shift_all>0)),numel(find(pos_shift_all~=0)),.5,'two');
pBinPosShiftAllNeg = myBinomTest(numel(find(pos_shift_all<0)),numel(find(pos_shift_all~=0)),.5,'two');

% hd
pBinHDShiftAllPosit = myBinomTest(numel(find(hd_shift_all>0)),numel(find(hd_shift_all~=0)),.5,'two');
pBinHDShiftAllNeg = myBinomTest(numel(find(hd_shift_all<0)),numel(find(hd_shift_all~=0)),.5,'two');

% spd
pBinSpdShiftAllPosit = myBinomTest(numel(find(spd_shift_all>0)),numel(find(spd_shift_all~=0)),.5,'two');
pBinSpdShiftAllNeg = myBinomTest(numel(find(spd_shift_all<0)),numel(find(spd_shift_all~=0)),.5,'two');


%% next, histograms for the three pooled pos hd and spd cells by shift

figure(4);
subplot(3,1,1);
histogram(pos_shift_all,'FaceColor','r','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Pooled Position Encoding Cells');
axis([ -50 50 0 length(pos_shift_all)]);
ax = gca;
ax.FontSize = 20;

subplot(3,1,2);
histogram(hd_shift_all,'FaceColor','c','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Pooled Head Direction Encoding Cells');
axis([ -50 50 0 length(hd_shift_all)]);
ax = gca;
ax.FontSize = 20;

subplot(3,1,3);
histogram(spd_shift_all,'FaceColor','b','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Pooled Speed Encoding Cells');
axis([ -50 50 0 length(spd_shift_all)]);
ax = gca;
ax.FontSize = 20;

%% plot PH model: position shifts versus head direction shifts

%first, add some jitter

b = -0.5; %upper limit for random numbers to create noise
a = 0.5; %lower limit for random numbers to create noise

noisePPosHd = (b-a).*rand(length(posHdMod),1) + a;
noiseHPosHd = (b-1).*rand(length(posHdMod),1) +a;

%add the jitter
posHdModJit = posHdMod;
posHdModJit(:,2) = posHdMod(:,2) + noisePPosHd;
posHdModJit(:,3) = posHdMod(:,3) + noiseHPosHd;

figure(5);
scatter(posHdModJit(:,2),posHdModJit(:,3),75,'filled')
xlabel('Shift in Position (PH Cells)','FontSize',20);
ylabel('Shift in Head Direction (PH Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% HS - same as above

%make noise
noiseHHdSpd = (b-a).*rand(length(hdSpdMod),1) + a;
noiseSHdSpd = (b-a).*rand(length(hdSpdMod),1) + a;

%add the jitter
hdSpdModJit = hdSpdMod;

hdSpdModJit(:,2) = hdSpdMod(:,2) + noiseHHdSpd;
hdSpdModJit(:,3) = hdSpdMod(:,3) + noiseSHdSpd;

figure(6);
scatter(hdSpdModJit(:,2),hdSpdModJit(:,3),75,'r','filled');
xlabel('Shift in Head Direction (HS Cells)','FontSize',20);
ylabel('Shift in Speed (HS Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% PS - same as above

%make noise
noisePPosSpd = (b-a).*rand(length(posSpdMod),1) + a;
noiseSPosSpd = (b-a).*rand(length(posSpdMod),1) + a;

% add the jitter
posSpdModJit = posSpdMod;

posSpdModJit(:,2) = posSpdMod(:,2) + noisePPosSpd;
posSpdModJit(:,3) = posSpdMod(:,3) + noiseSPosSpd;

%plot
figure(7);
scatter(posSpdModJit(:,2),posSpdModJit(:,3),75,'k','filled');
xlabel('Shift in Position (PS Cells)','FontSize',20);
ylabel('Shift in Speed (PS Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% PHS - need three separate plots
%plotting p versus h,  and h versus s, p versus s

% make some noise
noisePPosHdSpd = (b-a).*rand(length(posHdSpdMod),1) + a;
noiseHPosHdSpd = (b-a).*rand(length(posHdSpdMod),1) + a;
noiseSPosHdSpd = (b-a).*rand(length(posHdSpdMod),1) + a;

% add the jitters

posHdSpdModJit = posHdSpdMod;

posHdSpdModJit(:,2) = posHdSpdMod(:,2) + noisePPosHdSpd;
posHdSpdModJit(:,3) = posHdSpdMod(:,3) + noiseHPosHdSpd;
posHdSpdModJit(:,4) = posHdSpdMod(:,4) + noiseSPosHdSpd;

%plot these data!

figure(8);
scatter(posHdSpdModJit(:,2),posHdSpdModJit(:,3),75,'m','filled')
xlabel('Shift in Position (PHS Cells)','FontSize',20);
ylabel('Shift in Head Direction (PHS Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(9);
scatter(posHdSpdModJit(:,3),posHdSpdModJit(:,4),75,'k','filled');
xlabel('Shift in Head Direction (PHS Cells)','FontSize',20);
ylabel('Shift in Speed (PHS Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(10);
scatter(posHdSpdModJit(:,2),posHdSpdModJit(:,4),75,'c','filled')
xlabel('Shift in Position (PHS Cells)','FontSize',20);
ylabel('Shift in Speed (PHS Cells)','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% histograms for P only, H only, and Spd only!

figure(11);
subplot(3,1,1);
histogram(posMod(:,2),'FaceColor','g','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Position Only Encoding Cells')
axis([-50 50 0 length(posMod)])
ax = gca;
ax.FontSize = 20;

subplot(3,1,2);
histogram(hdMod(:,2),'FaceColor','b','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Head Direction Only Encoding Cells');
axis([-50 50 0 length(hdMod)])
ax = gca;
ax.FontSize = 20;

subplot(3,1,3);
histogram(spdMod(:,2),'FaceColor','r','BinWidth',2);
xlabel('Shift in Design Matrix','FontSize',20);
ylabel('Number of Cells','FontSize',20);
legend('Speed Only Encoding Cells');
axis([-50 50 0 length(spdMod)]);
ax = gca;
ax.FontSize = 20;

%% PCA on the total position cells

%first, position all

pos_tuning_all = [];
for k = encodingCells100'
    if ismember(1,variables_all{k}) && AllForwardFinal_Pval{k} < 0.05
        pos_tuning = tuning_curves_all{k}(1:101^2)*50; %multiplied by 50 to get into seconds, and not bins
        pos_tuning_all = [pos_tuning_all; pos_tuning'];
    end
end

% do PCA

PCA_matrixPos = pos_tuning_all;

% do pre-processing

% divide by range [personal pref.. divide std, nothing]
PCA_matrixPos = PCA_matrixPos./repmat(range(PCA_matrixPos,2),1,size(PCA_matrixPos,2));

% mean subtract
PCA_matrixPos = PCA_matrixPos - repmat(mean(PCA_matrixPos,2),1,size(PCA_matrixPos,2));
% do PCA
[pcsPos, proj_dataPos, eigenvalPos, ~, var_explained] = pca(PCA_matrixPos);

%make vector to show pos, neg, and zero values of the shift

%preallocate
shiftIDPos = zeros(length(pos_shift_all),1);

for i = 1:length(pos_shift_all)
    if pos_shift_all(i)>0
        shiftIDPos(i) = 1;
    elseif pos_shift_all(i)==0
        shiftIDPos(i) = 0;
    elseif pos_shift_all(i) < 0
        shiftIDPos(i) = -1;
    end
end


% plot projected data!
% projected data
figure(12)
scatter(proj_dataPos(:,1),proj_dataPos(:,2),75,shiftIDPos,'filled')
xlabel('PC 1','FontSize',20);
ylabel('PC 2','FontSize',20);
ax = gca;
ax.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTERING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%quantify clustering
%create the shift labels

shiftLabelsPos = zeros(length(pos_shift_all),1);
for i = 1:numel(pos_shift_all)
    if pos_shift_all(i)< 0
        shiftLabelsPos(i) = 1;
    else
        shiftLabelsPos(i) = 0;
    end
end

%call function to get knn

[distToNearNeighbPos,nullDistPos,pvalPosPCAClust] = findNearestNeighbor_shuffleLabels(proj_dataPos(:,1),proj_dataPos(:,2),shiftLabelsPos,5,1000);


%% Next, PCA on the HD cells

%get the tuning curves for head direction encoding cells
hd_tuning_all = [];

for k = encodingCells100'
    if ismember(2,variables_all{k}) && AllForwardFinal_Pval{k} <.05
        if ismember(1,variables_all{k})
            hd_tuning = tuning_curves_all{k}(10202:10301)*50;% multiplied by 50 to get these data into seconds
        else
            hd_tuning = tuning_curves_all{k}(1:100)*50;
        end
        hd_tuning_all(end+1,:) = hd_tuning';
    
    end
    
end

PCA_matrixHd = hd_tuning_all;

%preprocessing steps - divide by range

PCA_matrixHd = PCA_matrixHd./repmat(range(PCA_matrixHd,2),1,size(PCA_matrixHd,2));

%next, mean subtract

PCA_matrixHd = PCA_matrixHd - repmat(mean(PCA_matrixHd,2),1,size(PCA_matrixHd,2));

        
%do the PCA

[pcsHd, proj_dataHd, eigenvalHd, ~, ~] = pca(PCA_matrixHd);

%create a vector to ID the cells which shift in certain directions

%preallocate
shiftIDHD = zeros(length(hd_shift_all),1);

for i = 1:length(hd_shift_all)
    if hd_shift_all(i) > 0
        shiftIDHD(i) = 1;
    elseif hd_shift_all(i) == 0
        shiftIDHD(i) = 0;
    elseif hd_shift_all(i) < 0
        shiftIDHD(i) = -1;
    end
end

%plot the projected data
figure(13);
scatter(proj_dataHd(:,1),proj_dataHd(:,2),75,shiftIDHD,'filled');
xlabel('PC 1','FontSize',20);
ylabel('PC 2','FontSize',20);
ax = gca;
ax.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTERING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shiftLabelsHD = zeros(length(hd_shift_all),1);
for i = 1:numel(hd_shift_all)
    if hd_shift_all(i) < 0
        shiftLabelsHD(i) = 1;
    else
        shiftLabelsHD(i) = 0;
    end
end

%call the func

[distToNearNeighbHD,nullDistHD,pvalHDPCAClust] = findNearestNeighbor_shuffleLabels(proj_dataHd(:,1),proj_dataHd(:,2),shiftLabelsHD,5,1000);

%% last, do PCA on the spd cells

spd_tuning_all = [];

for k = encodingCells100'
    if ismember(3,variables_all{k}) && AllForwardFinal_Pval{k} < .05
        if numel(variables_all{k}) == 1 || (numel(variables_all{k}) == 2 && ismember(4,variables_all{k}))
            spd_tuning = tuning_curves_all{k}(1:60)*50; %if 3 only or 3 and 4
        elseif numel(variables_all{k}) == 2 && ismember(1,variables_all{k})
            spd_tuning = tuning_curves_all{k}(10202:10261)*50; %if its  1 3
        elseif numel(variables_all{k}) == 2 && ismember(2,variables_all{k})
            spd_tuning = tuning_curves_all{k}(101:160)*50; %if its 23
        elseif numel(variables_all{k}) == 3 && (ismember(1,variables_all{k}) && ismember(2,variables_all{k}))
            spd_tuning = tuning_curves_all{k}(10302:10361)*50; % if its 123
        elseif numel(variables_all{k}) == 3 && (ismember(1,variables_all{k}) && ismember(4,variables_all{k}))
            spd_tuning = tuning_curves_all{k}(10202:10261)*50; %if its 1 3 4
        elseif numel(variables_all{k}) == 3 && (ismember(2,variables_all{k}) && ismember(4,variables_all{k}))
            spd_tuning = tuning_curves_all{k}(101:160)*50; %if its 234
        elseif numel(variables_all{k}) == 4
            spd_tuning = tuning_curves_all{k}(10302:10361)*50; % if its 1234   
        end
        spd_tuning_all(end+1,:) = spd_tuning';
    end
end


%now, the PCA!

PCA_matrixSpd = spd_tuning_all;

%do the preprocessing - first, divide by the range

PCA_matrixSpd = PCA_matrixSpd./repmat(range(PCA_matrixSpd,2),1,size(PCA_matrixSpd,2));
    
%next, divide by the mean

PCA_matrixSpd = PCA_matrixSpd./repmat(mean(PCA_matrixSpd,2),1,size(PCA_matrixSpd,2));

%then, do the pca

[pcsSpd, proj_dataSpd, eigenvalSpd, ~, ~] = pca(PCA_matrixSpd);

%create vector identifying the cells which shift in certain directions

%preallocate
shiftIDSpd = zeros(length(spd_shift_all),1);

for i = 1:length(spd_shift_all)
    if spd_shift_all(i) > 0
        shiftIDSpd(i) = 1;
    elseif spd_shift_all(i) == 0
        shiftIDSpd(i) = 0;
    elseif spd_shift_all(i) < 0
        shiftIDSpd(i) = -1;
    end
end

%plot the projected data!

figure(14);
scatter(proj_dataSpd(:,1),proj_dataSpd(:,2),75,shiftIDSpd,'filled')
xlabel('PC 1','FontSize',20);
ylabel('PC 2','FontSize',20);
ax = gca;
ax.FontSize = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTERING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shiftLabelsSpd = zeros(length(spd_shift_all),1);
for i = 1:numel(spd_shift_all)
    if spd_shift_all(i)<0
        shiftLabelsSpd(i) = 1;
    else
        shiftLabelsSpd(i) = 0;
    end
end


%call the func
[distToNearNeighbSpd,nullDistSpd,pvalSpdPCAClust] = findNearestNeighbor_shuffleLabels(proj_dataSpd(:,1),proj_dataSpd(:,2),shiftLabelsSpd,5,1000);

%% Grid cell analysis - plot of fraction of grid and non-grid pos encoding cells which shift pos and neg

% Gridness score: Gridness Centre Removed
gridScore = xls_data(:,23);
GridCells100 = [];
for i = posModCellsAll
    if gridScore(i) > .4
        GridCells100(end+1) = i;
    end    
end


GridShift = zeros(length(GridCells100),1);

for i = 1:length(GridCells100)
    if cell2mat(AllPShifts{GridCells100(i)}(1)) < .05
        GridShift(i) = AllVarShifts{GridCells100(i)}(1);
    elseif cell2mat(AllPShifts{GridCells100(i)}(1)) > .05
        GridShift(i) = 0;
    end
end

noGridCells = setdiff(posModCellsAll,GridCells100);

noGridShift = zeros(length(noGridCells),1); %preallocate

for i = 1:length(noGridCells)
    if cell2mat(AllPShifts{noGridCells(i)}(1)) < .05
        noGridShift(i) = AllVarShifts{noGridCells(i)}(1);
    elseif cell2mat(AllPShifts{noGridCells(i)}(1)) > .05
        noGridShift(i) = 0;
    end
end

%given that a grid cell shifts, find the fraction of grid cells which shift
%pos and neg - do the same for non-grid pos encoding cells

%pos shift for grid
pGridPositShift = numel(find(GridShift>0))/numel(GridShift);

%neg shift for grid
pGridNegShift = numel(find(GridShift<0))/numel(GridShift);

%probaility of shifting given its a grid cell

pShiftGrid = numel(find(GridShift~=0))/numel(GridShift);

%same for non grid cells

%pos
pNoGridPositShift = numel(find(noGridShift>0))/numel(noGridShift);

%neg
pNoGridNegShift = numel(find(noGridShift<0))/numel(noGridShift);

%prob of shifting given not being a grid cell

pNoGridShift = numel(find(noGridShift~=0))/numel(noGridShift);

%plot these data

%grid data for shifting
grid_noGridData = ([pShiftGrid pGridPositShift pGridNegShift pNoGridShift pNoGridPositShift pNoGridNegShift]);

%labels
grid_noGridLabels = categorical({'Grid Cell Shift','Grid Cell Positive Shift','Grid Cell Negative Shift'...
                                 'Position-Encoding Non-Grid Cell Shift','Position-Encoding Non-Grid Cell Positive Shift'...
                                 'Position-Encoding Non-Grid Cell Negative Shift'});
grid_noGridLabels = reordercats(grid_noGridLabels,{'Grid Cell Shift','Grid Cell Positive Shift','Grid Cell Negative Shift'...
                                 'Position-Encoding Non-Grid Cell Shift','Position-Encoding Non-Grid Cell Positive Shift'...
                                 'Position-Encoding Non-Grid Cell Negative Shift'});

%plot these data

figure(15);
bar(grid_noGridLabels,grid_noGridData);
xlabel('Specific Category of Position Cell','FontSize',20);
ylabel('Probability','FontSize',20);
ax = gca;
ax.FontSize = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DIFFERENCES IN FRACTION OF SHIFTING BEHAVIOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%doing proportion test

%first, is the proportion of grid cells which neg shift sig larger than the
%proportion which shift positively

[zGridNegShift,pValGridNegShift] = prop_test(numel(find(GridShift<0)),numel(find(noGridShift<0)),numel(GridShift),numel(noGridShift));


% proportion test for positive shifts
[zGridPosShift,pValGridPosShift] = prop_test(numel(find(GridShift>0)),numel(find(noGridShift>0)),numel(GridShift),numel(noGridShift));


%% do the same procedure for border cells!

%find border cells

%loads the border score
borderScore = xls_data(:,45); 

%initialize
borderCells = [];

for i = encodingCells100'
    if borderScore(i) >.6
        borderCells(end+1) = i;
    end    
end

%transpose
borderCells = borderCells';

% find the non-border pos cells

noBorderCells = setdiff(posModCellsAll,borderCells);

%figure out shifting behavior

%preallocate
borderShift = zeros(length(borderCells),1);

for n = 1:numel(borderCells)
    if cell2mat(AllPShifts{borderCells(n)}(1)) < .05
        borderShift(n) = AllVarShifts{borderCells(n)}(1);  
    end
end

%no border shift

noBorderShift = zeros(length(noBorderCells),1);

for n = 1:numel(noBorderCells)
    if cell2mat(AllPShifts{noBorderCells(n)}(1))
        noBorderShift(n) = AllVarShifts{noBorderCells(n)}(1);
    end
end

%compute the probabilities

%prob of shift
pBorderShift = numel(find(borderShift~=0))/numel(borderShift);

%prob of pos shifft
pBorderPositShift = numel(find(borderShift>0))/numel(borderShift);

%prob of neg shift
pBorderNegShift = numel(find(borderShift<0))/numel(borderShift);

%same with no border
pNoBorderShift = numel(find(noBorderShift~=0))/numel(noBorderCells);

%pos shift
pNoBorderPositShift = numel(find(noBorderShift>0))/numel(noBorderShift);

%neg shift
pNoBorderNegShift = numel(find(noBorderShift<0))/numel(noBorderShift);

%plot these data

noBorder_borderData = [pBorderShift pBorderPositShift pBorderNegShift pNoBorderShift pNoBorderPositShift pNoBorderNegShift];

noBorder_borderLabels = categorical({'Border Cell Shift','Border Cell Positive Shift','Border Cell Negative Shift',...
                                      'Position-encoding non-Border Cell Shift','Position-encoding non-Border Cell Positive Shift','Position-encoding non-Border Cell Negative Shift'});
noBorder_borderLabels = reordercats(noBorder_borderLabels,{'Border Cell Shift','Border Cell Positive Shift','Border Cell Negative Shift',...
                                      'Position-encoding non-Border Cell Shift','Position-encoding non-Border Cell Positive Shift','Position-encoding non-Border Cell Negative Shift'});

%plot
figure(16);
bar(noBorder_borderLabels,noBorder_borderData)
xlabel('Category of Position Cell','FontSize',20);
ylabel('Probability','FontSize',20);
ax = gca;
ax.FontSize = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DIFFERENCES IN FRACTION OF SHIFTING BEHAVIOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first, look at negative shifting behavior
[zBorderNegShift,pValBorderNegShift] = prop_test(numel(find(borderShift<0)),numel(find(noBorderShift<0)),numel(borderShift),numel(noBorderShift));

%next look at positive shifting behavior
[zBorderPosShift,pValBorderPosShift] = prop_test(numel(find(borderShift>0)),numel(find(noBorderShift>0)),numel(borderShift),numel(noBorderShift));

%% do the same analysis for HD cells, looking at mean vector length

%HD score
meanVecLen = xls_data(:,43);

%initialize
hdScoreCells = [];

for i = hdModCellsAll
    if meanVecLen(i) > .2
        hdScoreCells(end+1) = i;      
    end
end

%transpose
hdScoreCells = hdScoreCells';

%find the hd shift

%preallocate
hdScoreShift = zeros(length(hdScoreCells),1);

for i = 1:numel(hdScoreCells)
    if ismember(1,variables_all{hdScoreCells(i)}) && ismember(2,variables_all{hdScoreCells(i)})
        if cell2mat(AllPShifts{hdScoreCells(i)}(2))< .05
            hdScoreShift(i) = AllVarShifts{hdScoreCells(i)}(2);
        end
    elseif ismember(2,variables_all{hdScoreCells(i)})
        if cell2mat(AllPShifts{hdScoreCells(i)}(1)) < .05
            hdScoreShift(i) = AllVarShifts{hdScoreCells(i)}(1);
        end
    else
        hdScoreShift(i) = 0;
    end
end

%find no score hd cells
hdNoScoreCells = setdiff(hdModCellsAll,hdScoreCells);

%find the shifts

%preallocate
hdNoScoreShift = zeros(numel(hdNoScoreCells),1);

for i = 1:numel(hdNoScoreCells)
    if ismember(1,variables_all{hdNoScoreCells(i)}) && ismember(2,variables_all{hdNoScoreCells(i)})
        if cell2mat(AllPShifts{hdNoScoreCells(i)}(2))< .05
            hdNoScoreShift(i) = AllVarShifts{hdNoScoreCells(i)}(2);
        end
    elseif ismember(2,variables_all{hdNoScoreCells(i)})
        if cell2mat(AllPShifts{hdNoScoreCells(i)}(1)) < .05
            hdNoScoreShift(i) = AllVarShifts{hdNoScoreCells(i)}(1);
        end
    else
        hdNoScoreShift(i) = 0;
    end
end

%compute the probabilities

pHdScoreShift = numel(find(hdScoreShift~=0))/numel(hdScoreShift);

pHdScorePositShift = numel(find(hdScoreShift>0))/numel(hdScoreShift);

pHdScoreNegShift = numel(find(hdScoreShift<0))/numel(hdScoreShift);

pHdNoScoreShift = numel(find(hdNoScoreShift~=0))/numel(hdScoreShift);

pHdNoScorePositShift = numel(find(hdNoScoreShift>0))/numel(hdNoScoreShift);

pHdNoScoreNegShift = numel(find(hdNoScoreShift<0))/numel(hdNoScoreShift);

%put the data together

hdScoreData = [pHdScoreShift pHdScorePositShift pHdScoreNegShift pHdNoScoreShift pHdNoScorePositShift pHdNoScoreNegShift];
hdScoreLabels = categorical({'HD Score Cell Shift','HD Score Cell Positive Shift','HD Score Cell Negative Shift', 'HD no Score Cell Shift','HD no Score Cell Positive Shift','HD no Score Cell Negative Shift'});
hdScoreLabels = reordercats(hdScoreLabels,{'HD Score Cell Shift','HD Score Cell Positive Shift','HD Score Cell Negative Shift','HD no Score Cell Shift','HD no Score Cell Positive Shift','HD no Score Cell Negative Shift'});

%plot data
figure(17);
bar(hdScoreLabels,hdScoreData);
xlabel('Category of Head Direction (By Score) Cell','FontSize',20);
ylabel('Probability','FontSize',20);
ax = gca;
ax.FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEED TO STATISTICALLY QUANTIFY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DIFFERENCES IN FRACTION OF SHIFTING BEHAVIOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use the binomial test to see if a HD cell is more likely to shift negatively
%or positively
[zHDScoreNegShift,pValHDScoreNegShift] = prop_test(numel(find(hdScoreShift<0)),numel(find(hdNoScoreShift<0)),numel(hdScoreShift),numel(hdNoScoreShift));
[zHDScorePositShift,pValHDScorePositShift] = prop_test(numel(find(hdScoreShift>0)),numel(find(hdNoScoreShift>0)),numel(hdScoreShift),numel(hdNoScoreShift));



%% exploratory analysis - looking at grid score and shift

figure(18);
scatter(gridScore(GridCells100),GridShift,75,'filled');
xlabel('Grid Score','FontSize',20);
ylabel('Shift in Design Matrix','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% conjunctive cells statisitcal indpendence analysis - first, PH

 %for PH cells: PH (pshift)/ph cells ; ph (h shift)/ph cells
% then doP(PH(pos)/p(PH(hshift))

PH_pPosShift = numel(find(posHdMod(:,2)~=0))/numel(posHdMod(:,2));

PH_pHDShift = numel(find(posHdMod(:,3)~=0))/numel(posHdMod(:,3));

%now, interested in the pPos shift given the pHD shift! compute the
%conditional prob

PH_pPos_gpHD = (PH_pPosShift*PH_pHDShift)/PH_pHDShift;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONSULT KIAH RE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% next, independence test for HS:

HS_pHDShift = numel(find(hdSpdMod(:,2)~=0))/numel(hdSpdMod(:,2));

HS_pSpdShift = numel(find(hdSpdMod(:,3)~=0))/numel(hdSpdMod(:,3));

%compute conditional probabililty

PH_pHD_gpSpd = (HS_pHDShift*HS_pSpdShift)/HS_pSpdShift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONSULT KIAH RE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Independence test for PS

PS_pPosShift = numel(find(posSpdMod(:,2)~=0))/numel(posSpdMod(:,2));

PS_pSpdShift = numel(find(posSpdMod(:,3)~=0))/numel(posSpdMod(:,3));

%compute conditional probability

PS_pPos_gpSpd = (PS_pPosShift*PS_pSpdShift)/PS_pSpdShift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONSULT KIAH RE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% exploratory analysis 2: look at running speed!

%get running speed score
spdScore = xls_data(:,end);

%first, get the spd cells running speed score

spd_score_spdModall = spdScore(spdModCellsAll);

%plot the scores with the shifts, in spd first

figure(19);
scatter(spd_score_spdModall,spd_shift_all,75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in All Running Speed Encoding Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% look at how running speed affects speed only

spd_score_spdMod = spdScore(spdMod(:,1));

figure(20);
scatter(spd_score_spdMod,spdMod(:,2),75,'filled')
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Running Speed Only Encoding Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% look at how running speed modulates postion - so, PS!

%get scores for these cells
spd_score_posSpd = spdScore(posSpdMod(:,1));

figure(21);
scatter(spd_score_posSpd,posSpdMod(:,2),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Position in PS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(22);
scatter(spd_score_posSpd,posSpdMod(:,3),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Running Speed in PS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% next, HS cells

figure(23);
hdSpdModSpdScore = spdScore(hdSpdMod(:,1));
scatter(hdSpdModSpdScore,hdSpdMod(:,2),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Head Direction in HS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(24);
scatter(hdSpdModSpdScore,hdSpdMod(:,3),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Running Speed in HS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% last, PHS!

posHdSpdModSpdScore = spdScore(posHdSpdMod(:,1));

figure(25);
scatter(posHdSpdModSpdScore,posHdSpdMod(:,2),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Position in PHS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(25);
scatter(posHdSpdModSpdScore,posHdSpdMod(:,3),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Head Direction in PHS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

figure(26);
scatter(posHdSpdModSpdScore,posHdSpdMod(:,4),75,'filled');
xlabel('Speed Score','FontSize',20);
ylabel('Shift in Running Speed in PHS Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% for grid, border, hd analysis - find numbers in each

figure(27);
numGridBorderHDSc = [numel(GridCells100) numel(noGridCells) numel(borderCells) numel(noBorderCells) numel(hdScoreCells) numel(hdNoScoreCells)];

labelsGridBorderHDSc = categorical({'Grid Cells','Position Encoding non-Grid Cells','Border Cells','Position Encoding non-Border Cells','Head Direction Cells (Score > .2)','Head Direction Cells (Score < .2)'});

labelsGridBorderHDSc = reordercats(labelsGridBorderHDSc,{'Grid Cells','Position Encoding non-Grid Cells','Border Cells','Position Encoding non-Border Cells','Head Direction Cells (Score > .2)','Head Direction Cells (Score < .2)'});

bar(labelsGridBorderHDSc,numGridBorderHDSc);
xlabel('Cell Type','FontSize',20);
ylabel('Number of Cells','FontSize',20);
ax = gca;
ax.FontSize = 20;

%% more analysis - first, look at model fit as a function of shifts
%interested specifically in how the single variables within 








%tuning curves for position shifts

%speed tuning curves for positive and negative shifts

% one variable shifting for conjunctive cell..explore how the model fits
% for each variable may differ and if this is related to the interesting
% shifting behavior

%first, the PH cells

%initialize matrices for finding the model order (i.e., the order in which
%the variables were added to the model
%reset the name of these
ph = [1 2];
allFits = AllvarAllmodelFits;

%preallocate
orderPosHdMod = zeros(length(posHdMod),2);


for i = 1:numel(posHdMod(:,1))
    for j = 1:numel(ph)
        [~,orderPosHdMod(i,j)] = max(allFits{posHdMod(i,1)}(1:3,j));
    end
end


% next,try to see if there is relationship between the order of addition
% to the model and the shifts

%first position - plot the order in which it was added and the observed
%shift

figure(28);
subplot(2,1,1);
histogram(posHdMod(find(orderPosHdMod(:,1)==1),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdMod)])
title('Pos = 1, HD = 2')
subplot(2,1,2);
histogram(posHdMod(find(orderPosHdMod(:,1)==1),3),'BinWidth',1.5,'FaceColor','k');
axis([ -50 50 0 length(posHdMod)])
title('Pos = 1, HD = 2')

figure(29);
subplot(2,1,1);
histogram(posHdMod(find(orderPosHdMod(:,1)==2),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdMod)])
title('Pos = 2, HD = 1')
subplot(2,1,2);
histogram(posHdMod(find(orderPosHdMod(:,1)==2),3),'BinWidth',1.5,'FaceColor','k');
axis([ -50 50 0 length(posHdMod)])
title('Pos = 2, HD = 1')

%% same thing for the PS model

allFits = AllvarAllmodelFits;

ps = [1 3];

orderPosSpdMod = zeros(length(posSpdMod),2);

for i = 1:numel(posSpdMod(:,1))
    for j = 1:numel(ps)
        [~,orderPosSpdMod(i,j)] = max(allFits{posSpdMod(i,1)}(1:3,j));
    end
end

% plot to see relationship

figure(30);
subplot(2,1,1);
histogram(posSpdMod(find(orderPosSpdMod(:,1)==1),2),'BinWidth',1.5);
axis([ -50 50 0 length(posSpdMod)])
title('Pos = 1, Spd = 2');
subplot(2,1,2);
histogram(posSpdMod(find(orderPosSpdMod(:,1)==1),3),'BinWidth',1.5,'FaceColor','c');
axis([ -50 50 0 length(posSpdMod)])
title('Pos = 1, Spd = 2');

figure(31);
subplot(2,1,1);
histogram(posSpdMod(find(orderPosSpdMod(:,1)==3),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdMod)])
title('Pos = 2, Spd = 1');
subplot(2,1,2);
histogram(posSpdMod(find(orderPosSpdMod(:,1)==3),3),'BinWidth',1.5,'FaceColor','g');
axis([ -50 50 0 length(posHdMod)])
title('Pos = 2, Spd = 1');

%% same but for HS!!!!!!!!

allFits = AllvarAllmodelFits;
hs = [2 3];

orderHdSpdMod = zeros(length(hdSpdMod),2);

for i = 1:numel(hdSpdMod(:,1))
    for j = 1:numel(hs)
        [~,orderHdSpdMod(i,j)] = max(allFits{hdSpdMod(i,1)}(1:3,j));
    end
end

%plot to visualize

figure(32);
subplot(2,1,1);
histogram(hdSpdMod(find(orderHdSpdMod(:,1)==2),2),'BinWidth',1.5)
axis([ -50 50 0 length(posHdMod)])
title('HD = 1, Spd = 2');
subplot(2,1,2);
histogram(hdSpdMod(find(orderHdSpdMod(:,1)==2),3),'BinWidth',1.5,'FaceColor','r');
axis([ -50 50 0 length(posHdMod)])
title('HD = 1, Spd = 2');


figure(33);
subplot(2,1,1);
histogram(hdSpdMod(find(orderHdSpdMod(:,1)==3),2),'BinWidth',1.5);
title('HD = 2, Spd = 1');
axis([ -50 50 0 length(posHdMod)])
subplot(2,1,2);
histogram(hdSpdMod(find(orderHdSpdMod(:,1)==3),3),'BinWidth',1.5,'FaceColor','k');
axis([ -50 50 0 length(posHdMod)])
title('HD = 1, Spd = 2');

%% same thing for the PHS

allFits = AllvarAllmodelFits;

phs = [1 2 3];

orderPosHdSpdMod = zeros(length(posHdSpdMod),3);

for i = 1:numel(posHdSpdMod(:,1))
    for j = 1:numel(phs)
        [~,orderPosHdSpdMod(i,j)] = max(allFits{posHdSpdMod(i,1)}(1:3,j));
    end
end

% plot to visualize

%first compare the position to head direction
figure(34);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==1),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 1, HD = 2');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==1),3),'BinWidth',1.5,'FaceColor','b');
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 1, HD = 2');

figure(35);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==2),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 2, HD = 1');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==2),3),'BinWidth',1.5,'FaceColor','b');
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 2, HD = 1');

%next, compare position and spd

figure(36);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==1),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 1, Spd = 2');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==1),4),'BinWidth',1.5,'FaceColor','c');
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 1, Spd = 2');

figure(37);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==3),2),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 2, Spd = 1');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==3),4),'BinWidth',1.5,'FaceColor','c');
axis([ -50 50 0 length(posHdSpdMod)])
title('Pos = 2, Spd = 1');

% last, look at HS!

figure(38);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1) == 2),3),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('HD = 1, Spd = 2');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1) == 2),4),'BinWidth',1.5,'FaceColor','c');
axis([ -50 50 0 length(posHdSpdMod)])
title('HD = 1, Spd = 2');

figure(39);
subplot(2,1,1);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1) == 3),3),'BinWidth',1.5);
axis([ -50 50 0 length(posHdSpdMod)])
title('HD = 2, Spd = 1');
subplot(2,1,2);
histogram(posHdSpdMod(find(orderPosHdSpdMod(:,1)==3),4),'BinWidth',1.5,'FaceColor','g');
axis([ -50 50 0 length(posHdSpdMod)])
title('HD = 2, Spd = 1');

%% next, look at tuning curves for the position cells which shift versus those which don't shift

%more specifically, we are interested in non-grid cells which shifted
%versus did not shift:

%first, the non-grid cells which do shift
%counter
c = 0;
for i = 1:numel(noGridCells)
    if noGridShift(i) ~=0
       c = c+1; %advance counter
       h = figure;
       imagesc(reshape(pos_tuning_all(i,:),101,101));
       saveas(h,sprintf('FIG_noGridShift%d.png',c));
    end
end

%do the same thing for the cells which don't shift
%reset the counter
c = 0;
for i = 1:numel(noGridCells)
    if noGridShift(i) == 0
        c = c+1; %advance the counter
        h = figure;
        imagesc(reshape(pos_tuning_all(i,:),101,101))
        saveas(h,sprintf('FIG_noGridNoShift%d.png',c));
    end
end

%% last, look at tuning curves for speed cells which shifted positive versus negative

%first shifting speed cells

%reset the counter
for i = 1:numel(spdModCellsAll)
    if spd_shift_all(i) ~=0
        c = c+1; %advance the counter
        h = figure;
        plot(spd_tuning_all(i,:),'k','linewidth',2);
        saveas(h,sprintf('FIG_SpdShift%d.png',c));    
    end
end

%same thing for non-shifting speed cells

%reset the count 
c =0;
for i = 1:numel(spdModCellsAll)
    if spd_shift_all(i)==0
        c = c+1;
        h = figure;
        plot(spd_tuning_all(i,:),'k','LineWidth',2);
        saveas(h,sprintf('FIG_SpdNoShift%d.png',c));
    end    
end


%% make matrix for kiah!!!

%has all of the cells and their shifts for pos, hd, and spd - here we will
%get rid of any shifts which weren't significant by setting them to zero

%preallocate for new shifts - this is for significant ones only
SigAllVarShifts = cell(size(AllVarShifts));

%preallocate for the new matrix
allCellsModsShifts = zeros(numel(encodingCells),3);

%fill up the matrix with all of the shifts

for i = encodingCells'
    %first the position encoding
    if ismember(1,variables_all{i}) && AllForwardFinal_Pval{i} < .05 && cell2mat(AllPShifts{i}(find(variables_all{i}==1))) < .05
        allCellsModsShifts(i,1) = AllVarShifts{i}(find(variables_all{i} == 1));
    else
        allCellsModsShifts(i,1) = 0;
    end
    %next hd encoding
    if ismember(2,variables_all{i}) && AllForwardFinal_Pval{i} < .05 && cell2mat(AllPShifts{i}(find(variables_all{i}==2))) < .05
        allCellsModsShifts(i,2) = AllVarShifts{i}(find(variables_all{i} == 2));
    else
        allCellsModsShifts(i,2) = 0;
    end
    %last the spd
    if ismember(3,variables_all{i}) && AllForwardFinal_Pval{i} < .05 && cell2mat(AllPShifts{i}(find(variables_all{i}==3))) < .05
        allCellsModsShifts(i,3) = AllVarShifts{i}(find(variables_all{i} == 3));
    else
        allCellsModsShifts(i,3) = 0;
    end
end













 








