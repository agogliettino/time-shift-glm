clear all; close all; clc %clear the workspace


%fig 1 cell type verusus number of cells
%fig 2 - cell type versus number of cells - divide bars into positive
%shifts, neg shifts, or zero shifts
%fig 3 - then do the all-position, all hd, all spd verusus # cells (divide
%up bars)
%fig 4 - do position all, hd all, spd all versus cell number, dividing up
%the bar by type of shift
%for all them - figure out if the distribution of shifts for each type is
%diff than chance (binomial)
%fig 5 - 8: shifts for pos all versus cells; hd versus cells; spd versus cells
%fig 9 - PH cells plot pos shift versus HD shift
%fig 10 - HS cells: hd versus spd 
%fig 11 - PS (pos vs spd)
%fig 12 - PHS: plot ph, hs, ps shifts (3 plots)
%fig 13-15 - shift histograms for p only, h only, and s only
%fig 16 - 18: PCA plots for pos, hd, spd (make the shift type categorical)
%and p val for pos and neg shift clustering (two p vals per plot bc of pos
%and neg shifts)
%fig 19 - grid score: plot of fraction of cells that shift pos and neg for
%both grid cells and non grid cells (pos cells) - fig 20 - compute p val for
%different fractions of posit and neg shifts for grid cells vers non grid
%cells
%fig 20, 21 border cells (above 0.6); head direction: mean vector length
%(.02) do same analysis as for grid cells

%exploratory - avg grid score and shift - probability
%conj cells  - for PH cells: PH (pshift)/ph cells ; ph (h shift)/ph cells
% then doP(PH(pos)/p(PH(hshift))

%do the same for the other models, too - wanna find out if they're
%independent 

% for the model:
%step 1 - find the shift
%step 2 - fit and compare
%assump - variables are indpenedent and the param wont change too much


%% load data

%load the cells of interest
load('encodingCells.mat');

load('do_glm_output_112817.mat'); %loads the data into the workspace - this contains all of the cells


xls_data = xlsread('Wildtype_Mice_Database_v5.xlsx');
boxSizeCell = xls_data(:,13);
cellboxSize100 = find(boxSizeCell == 100);
encodingCells100 = intersect(encodingCells,cellboxSize100);

%% find the pos, HD, and spd cells
%these are cells which encode these variables
%first, the HD
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

%% total number of shifting cells

%number of position shift
numPosShift = find(numel(pos_shift_all ~=0));

%specifically, look at grid cells!

% Gridness score: Gridness Centre Removed
gridScore = xls_data(:,23);
GridCells = find(gridScore > 0.4);

%get grid cells with box == 100
GridCells100 = intersect(GridCells,encodingCells100);

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

%figure out the probability of shifitng for grid cells out of all of the
%pos cells which shift

pGridShift = numel(find(GridShift~=0))/numel(GridCells100);

%tuning curves for non grid cells and statistically quantify



%probability of shifting 

pNoGridShift = numel(find(noGridShift~=0))/numel(noGridCells);

%hd
numHdShift = numel(find(hd_shift_all ~=0));

%spd
numSpdShift = numel(find(spd_shift_all ~=0));

% plot the shifts for all the different cell types

figure(1);
bar(posModCellsAll,pos_shift_all,'FaceColor','k');
hold on
bar(hdModCellsAll,hd_shift_all,'FaceColor','m');
hold on
bar(spdModCellsAll,spd_shift_all,'FaceColor','c');
legend('Pos','HD','Spd');
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%% quantify tendency of a variable to shift for a cell which encodes it

%for position shift
pPosShift = numel(find(pos_shift_all < 0))/numel(find(pos_shift_all>=0));

%for hd shift

pHdShift = numel(find(hd_shift_all < 0))/numel(find(hd_shift_all >=0));

%for the spd shift

pSpdShift = numel(find(spd_shift_all < 0))/numel(find(spd_shift_all>=0));

%% separate cells which are specifically model encoding

%first, group cells based on the different possible models

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


%% get rid of the empty cells

%convert from cell array to matrix

posMod = cell2mat(posMod);
hdMod = cell2mat(hdMod);
spdMod = cell2mat(spdMod);
posHdMod = cell2mat(posHdMod);
hdSpdMod = cell2mat(hdSpdMod);
posSpdMod = cell2mat(posSpdMod);
posHdSpdMod = cell2mat(posHdSpdMod);

%plot these data

%first the position models
figure(2);
histogram(posMod(:,2),'FaceColor','r','BinWidth',1)
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Pos Mod');

%hd
figure(3);
histogram(hdMod(:,2),'FaceColor','g','BinWidth',1)
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Hd Mod');

%spd - had nothing interesting ie. all zero shifts

%next, for the cells encoding multiple variables - first, the pos HD model

figure(4);
subplot(2,1,1)
histogram(posHdMod(:,2),'FaceColor','b','BinWidth',1); %pos 
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Pos Shifts');
axis([-50 50 0 length(posHdMod(:,2))])
subplot(2,1,2)
histogram(posHdMod(:,3),'FaceColor','r','BinWidth',1); %HD
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('HD Shifts');
axis([-50 50 0 length(posHdMod(:,3))])

%next, the hdSpd
figure(5);
subplot(2,1,1)
histogram(hdSpdMod(:,2),'FaceColor','g','BinWidth',1); %hd 
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('HD Shift');
axis([-50 50 0 length(hdSpdMod(:,2))])
subplot(2,1,2)
histogram(hdSpdMod(:,3),'FaceColor','b','BinWidth',1); %spd
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Spd Shift');
axis([-50 50 0 length(hdSpdMod(:,3))])

%posSpd
figure(6);
subplot(2,1,1)
histogram(posSpdMod(:,2),'FaceColor','r','BinWidth',1); %pos 
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Pos Shift');
axis([-50 50 0 length(posSpdMod(:,2))])
subplot(2,1,2)
histogram(posSpdMod(:,3),'FaceColor','g','BinWidth',1); %Spd
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Spd Shift');
axis([-50 50 0 length(posSpdMod(:,3))])

%full model
figure(7);
subplot(3,1,1)
histogram(posHdSpdMod(:,2),'FaceColor','c','BinWidth',1); % position
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Pos Shift');
axis([-50 50 0 length(posHdSpdMod(:,2))])
subplot(3,1,2)
histogram(posHdSpdMod(:,3),'FaceColor','m','BinWidth',1); %HD
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('HD Shift');
axis([-50 50 0 length(posHdSpdMod(:,3))])
subplot(3,1,3);
histogram(posHdSpdMod(:,4),'FaceColor','g','BinWidth',1); %spd
xlabel('Shift in Design Matrix');
ylabel('Number of cells');
legend('Spd Shift');
axis([-50 50 0 length(posHdSpdMod(:,4))])
keyboard

%% probability of shifting given an encoded variable
% prob of shifting given a cell 

%position first
pPosShift = numel(find(pos_shift_all ~=0))/numel(pos_shift_all);

%hd
pHdShift = numel(find(hd_shift_all ~=0))/numel(hd_shift_all);

%spd
pSpdShift = numel(find(spd_shift_all ~= 0))/numel(spd_shift_all);

%% probability of shifting in a cerain direction, given an encoded variable

%position
pPosShiftPos = numel(find(pos_shift_all > 0))/numel(find(pos_shift_all ~=0));
pPosShiftNeg = numel(find(pos_shift_all < 0))/numel(find(pos_shift_all ~=0));

%hd
pHdShiftPos = numel(find(hd_shift_all > 0))/numel(find(hd_shift_all ~=0));
pHdShiftNeg = numel(find(hd_shift_all < 0))/numel(find(hd_shift_all ~=0));

%spd
pSpdShiftPos = numel(find(spd_shift_all > 0))/numel(find(spd_shift_all ~=0));
pSpdShiftNeg = numel(find(spd_shift_all < 0))/numel(find(spd_shift_all ~=0));

%% PCA 

run doPCA_kh.m

%% plot tuning curves for grid cells which shift versus those which dont

%define grid cells which do shift
GridCellsShift = GridCells100(find(GridShift~=0));

[~,GridCellsShiftIdx] = intersect(posModCellsAll,GridCellsShift);

%define grid cells which do not shift

GridCellsNoShift = GridCells100(find(GridShift == 0));

[~,GridCellsNoShiftIdx] = intersect(posModCellsAll,GridCellsNoShift);

[~,noGridCellsIdx] = setdiff(posModCellsAll,noGridCells);

%first, tuning curves for the grid cells which shift
c = 1;
for n = GridCellsShiftIdx'
    figure(c);
    imagesc(reshape(pos_tuning_all(n,:),101,101))
    c = c+1; %advance the counter
end

%next, tuning curves for the grid cells which don't shift
for n = GridCellsNoShiftIdx'
    figure(c);
    imagesc(reshape(pos_tuning_all(n,:),101,101))
    c = c + 1;
end

%look at tuning curves for non grid cells - sanity check :)

for n = noGridCellsIdx'
    figure(c);
    imagesc(reshape(pos_tuning_all(n,:),101,101));
    c = c+1;
end





















