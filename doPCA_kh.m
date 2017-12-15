% doPCA

%% clear the workspace
clear all; close all; clc

%loads the cells
load('encodingCells.mat'); 

%load the data on which the glm was run
load('do_glm_output_112817.mat');

xls_data = xlsread('Wildtype_Mice_Database_v5.xlsx');
boxSizeCell = xls_data(:,13);
cellboxSize100 = find(boxSizeCell == 100);
encodingCells100 = intersect(encodingCells,cellboxSize100);

% Gridness score: Gridness Centre Removed
gridScore = xls_data(:,23);
GridCells = find(gridScore > 0.4);

%get grid cells with box == 100
GridCells100 = intersect(GridCells,encodingCells100);

%next compute if this is significant
%plot this and non grid
%look at tuning curves of non-grid posiotn cells whcih exhibit position
%shifts

%% collect all the tuning curves for position encoding cells

pos_tuning_all = [];
pos_shift_all = [];
posModCellsAll = [];
for k = encodingCells100'
    if ismember(1,variables_all{k}) && AllForwardFinal_Pval{k} < 0.05
        posModCellsAll(end+1) = k; 
        pos_tuning = tuning_curves_all{k}(1:101^2)*50; %multiplied by 50 to get into seconds, and not bins
        pos_tuning_all = [pos_tuning_all; pos_tuning'];
        
        pval = AllPShifts{k}{1};
        if pval < 0.05
            pos_shift_all = [pos_shift_all; AllVarShifts{k}(1)];
        else
            pos_shift_all = [pos_shift_all; 0];
        end
    end
end
%% do PCA

PCA_matrixPos = pos_tuning_all;

% do pre-processing

% divide by range [personal pref.. divide std, nothing]
PCA_matrixPos = PCA_matrixPos./repmat(range(PCA_matrixPos,2),1,size(PCA_matrixPos,2));

% mean subtract
PCA_matrixPos = PCA_matrixPos - repmat(mean(PCA_matrixPos,2),1,size(PCA_matrixPos,2));
% do PCA
[pcs, proj_data, eigenval, ~, var_explained] = pca(PCA_matrixPos);

% look at eigenvalues
figure(1)
stem(eigenval)

% look at the eigenvectors / pcs --- not interesting generally, but useful
% for a sanity check
figure(2)
subplot(1,3,1)
imagesc(reshape(pcs(:,1),101,101))
subplot(1,3,2)
imagesc(reshape(pcs(:,2),101,101))
subplot(1,3,3)
imagesc(reshape(pcs(:,3),101,101))

% projected data
figure(3)
scatter(proj_data(:,1),proj_data(:,2),75,pos_shift_all,'filled')

%% test if clustering between cells which shifted negatively vs cells which did not shift/shift positively was significant

%first compute the relevant labels 

%label to test between negative vs pos cells/non shift cells


shiftingCells = [];
noShiftingCells = [];

for i = 1:length(posModCellsAll)
    if pos_shift_all(i) < 0
        shiftingCells(end+1) = posModCellsAll(i);
    elseif pos_shift_all(i) >= 0 
        noShiftingCells(end+1) = posModCellsAll(i);
    end
end

%create a label for assignment scrambling

%preallocate
shiftLabels = zeros(length(posModCellsAll),1);

for i = 1:numel(posModCellsAll)
    if ismember(posModCellsAll(i),shiftingCells)
        shiftLabels(i) = 1;
    elseif ismember(posModCellsAll(i),noShiftingCells)
        shiftLabels(i) = 0;
    end
end

% call function to determine if clustering is significant (1) or not (0)

[distToNearNeighb,nullDist,pval] = findNearestNeighbor_shuffleLabels(proj_data(:,1),proj_data(:,2),shiftLabels,5,1000);



%% test grid cell versus non-grid position cell shifting!

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
keyboard
%figure out the probability of shifitng for grid cells out of all of the
%pos cells which shift

pGridShift = numel(find(GridShift~=0))/numel(GridCells100);

%probability of shifting 

pNoGridShift = numel(find(noGridShift~=0))/numel(noGridCells);


%now, do pca for hd 

%get tuning curves for all significant hd cells

hd_tuning_all = [];
hd_shift_all = [];
hdModCellsAll = [];

for k = encodingCells100'
    if ismember(2,variables_all{k}) && AllForwardFinal_Pval{k} <.05
        hdModCellsAll(end+1) = k;
        if ismember(1,variables_all{k})
            hd_tuning = tuning_curves_all{k}(10202:10301)*50;% multiplied by 50 to get these data into seconds
        else
            hd_tuning = tuning_curves_all{k}(1:100);
        end
        hd_tuning_all(end+1,:) = hd_tuning';
   
    
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
        elseif numel(variables_all{k}) == 2 && (isequal(variables_all{k},[2 3]) || isequal(variables_all{k},[2 4]))
            pval = AllPShifts{k}{1};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(1)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 3 && (isequal(variables_all{k},[1 2 3]) || isequal(variables_all{k},[1 2 4]))
            pval = AllPShifts{k}{2};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(2)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 3 && isequal(variables_all{k},[2 3 4])
            pval = AllPShifts{k}{1};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(1)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        elseif numel(variables_all{k}) == 4
            pval = AllPShifts{k}{2};
            if pval < 0.05
                hd_shift_all = [hd_shift_all; AllVarShifts{k}(2)];
            else
                hd_shift_all = [hd_shift_all; 0];
            end
        end
    
    end
    
end
keyboard
%do PCA on the hd model cells

PCA_matrixHd = hd_tuning_all;

%preprocessing steps - divide by range

PCA_matrixHd = PCA_matrixHd./repmat(range(PCA_matrixHd,2),1,size(PCA_matrixHd,2));

%next, mean subtract

PCA_matrixHd = PCA_matrixHd - repmat(mean(PCA_matrixHd,2),1,size(PCA_matrixHd,2));

        
%do the PCA

[pcsHd, proj_dataHd, eigenvalHd, ~, ~] = pca(PCA_matrixHd);

%eigenvalues

figure(4);
stem(eigenval)

%plot eigenvectors and Pcs

figure(5);
subplot(1,3,1);
plot(pcsHd(:,1),'k','linewidth',2);
subplot(1,3,2);
plot(pcsHd(:,2),'k','linewidth',2);
subplot(1,3,3);
plot(pcsHd(:,3),'k','linewidth',2);

%plot the projected data

figure(6);
scatter(proj_dataHd(:,1),proj_dataHd(:,2),75,hd_shift_all,'filled');

%check for a significance of clusters

shiftingCells = [];
noShiftingCells = [];

for i = 1:length(hdModCellsAll)
    if hd_shift_all(i) < 0
        shiftingCells(end+1) = hdModCellsAll(i);
    elseif hd_shift_all(i) >= 0 
        noShiftingCells(end+1) = hdModCellsAll(i);
    end
end

%create a label for assignment scrambling

%preallocate
shiftLabels = zeros(length(hdModCellsAll),1);

for i = 1:numel(hdModCellsAll)
    if ismember(hdModCellsAll(i),shiftingCells)
        shiftLabels(i) = 1;
    elseif ismember(hdModCellsAll(i),noShiftingCells)
        shiftLabels(i) = 0;
    end
end

%call the function to get sig

[distToNearNeighb,nullDist,pval] = findNearestNeighbor_shuffleLabels(proj_dataHd(:,1),proj_dataHd(:,2),shiftLabels,5,1000);


%last, the speed analysis
spd_tuning_all = [];
spd_shift_all = [];
spdModCellsAll = [];

%need to consider different conditions to evaluate if the cell has speed
%(and significantly) in its model

for k = encodingCells100'
    if ismember(3,variables_all{k}) && AllForwardFinal_Pval{k} < .05
        spdModCellsAll(end+1) = k;
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

%now, the PCA!

PCA_matrixSpd = spd_tuning_all;

%do the preprocessing - first, divide by the range

PCA_matrixSpd = PCA_matrixSpd./repmat(range(PCA_matrixSpd,2),1,size(PCA_matrixSpd,2));
    
%next, divide by the mean

PCA_matrixSpd = PCA_matrixSpd./repmat(mean(PCA_matrixSpd,2),1,size(PCA_matrixSpd,2));

%then, do the pca

[pcsSpd, proj_dataSpd, eigenvalSpd, ~, ~] = pca(PCA_matrixSpd);

%plot the eigen values

figure(7);
stem(eigenvalSpd)

%plot the pcs/eigenvectors

figure(8);
subplot(1,3,1);
plot(pcsSpd(:,1),'k','linewidth',2);
subplot(1,3,2);
plot(pcsSpd(:,2),'k','linewidth',2);
subplot(1,3,3);
plot(pcsSpd(:,3),'k','linewidth',2);

%plot the projected data

figure(9);
scatter(proj_dataSpd(:,1),proj_dataSpd(:,2),75,spd_shift_all,'filled')


%do analysis to figure out if clustering is sig


shiftingCells = [];
noShiftingCells = [];

for i = 1:length(spdModCellsAll)
    if spd_shift_all(i) < 0
        shiftingCells(end+1) = spdModCellsAll(i);
    elseif spd_shift_all(i) >= 0 
        noShiftingCells(end+1) = spdModCellsAll(i);
    end
end

%create a label for assignment scrambling

%preallocate
shiftLabels = zeros(length(hdModCellsAll),1);

for i = 1:numel(hdModCellsAll)
    if ismember(hdModCellsAll(i),shiftingCells)
        shiftLabels(i) = 1;
    elseif ismember(hdModCellsAll(i),noShiftingCells)
        shiftLabels(i) = 0;
    end
end

%call the function to get sig

[distToNearNeighb,nullDist,pval] = findNearestNeighbor_shuffleLabels(proj_dataHd(:,1),proj_dataHd(:,2),shiftLabels,5,1000);

    





