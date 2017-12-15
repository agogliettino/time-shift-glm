function [meanEdShift,meanEdNoShift,meanEdShiftScr,meanEdNoShiftScr,shiftLabels] = nearestNeighborwShuff(modShiftAll,cellsModelAll,dataRedProj)
% to test the effect's significance
% divide points into 2 bins: those that don't shift (or positively; >=0), and ones that shift
% negatively (< 0)
% for every point: compute avg distance to N nearest neighbors (5?)
% scramble assigments. recompute avg distance to N nearest neighbors
% compare data to shuffled distribution of avg distances
%first, for position cells - divide cells up into shifting versus
%non-shifting cells

cellsModelAll = cellsModelAll'; %transpose this

shiftingCells = [];
noShiftingCells = [];

for i = 1:length(cellsModelAll)
    if modShiftAll(i) < 0
        shiftingCells(end+1) = cellsModelAll(i);
    elseif modShiftAll(i) >= 0 
        noShiftingCells(end+1) = cellsModelAll(i);
    end
end

%create a label for assignment scrambling

%preallocate
shiftLabels = zeros(length(cellsModelAll),1);

for i = 1:numel(cellsModelAll)
    if ismember(cellsModelAll(i),shiftingCells)
        shiftLabels(i) = 1;
    elseif ismember(cellsModelAll(i),noShiftingCells)
        shiftLabels(i) = 0;
    end
end


%show where there are ones

shiftingIdx = find(shiftLabels == 1);
noShiftingIdx = find(shiftLabels == 0);
%preallocate for the means


%create matrix for projected data - ie the x and ys

proj_dataMat = [dataRedProj(:,1) dataRedProj(:,2)];

%for all the cells, find the 5 nearest neighbors, and their distances!

%preallocate for  both the neighbors and distances
nbShift = cell(numel(shiftingIdx),1); % neighboring cells
edShift = cell(numel(shiftingIdx),1); %euclidean distance

%do it first for the shifting
for n = 1:numel(shiftingIdx)
    [nbShift{n},edShift{n}] = knnsearch(proj_dataMat(shiftingIdx,:),proj_dataMat(shiftingIdx(n),:),'k',6);
end

%preallocate for  both the neighbors and distances
nbNoShift = cell(numel(noShiftingIdx),1); % neighboring cells
edNoShift = cell(numel(noShiftingIdx),1); %euclidean distance
%next, do it for the no shifting

for n = 1:numel(noShiftingIdx)
    [nbNoShift{n},edNoShift{n}] = knnsearch(proj_dataMat(noShiftingIdx,:),proj_dataMat(noShiftingIdx(n),:),'k',6);   
end

%compute the mean for both now 

meanEdShift = zeros(numel(shiftingIdx),1); %preallocate

for m = 1:length(edShift)
    meanEdShift(m) = mean(edShift{m}(2:6)); %dont want to include the self-distance
end

meanEdNoShift = zeros(numel(noShiftingIdx),1); %preallocate

for m = 1:length(edNoShift)
    meanEdNoShift(m) = mean(edNoShift{m}(2:6));
end

%now, shuffle the assignments, and look at how this distance changes
shiftLabelsScr = shiftLabels(randperm(numel(shiftLabels)));

%with this, recompute the shifting and no shfitingIDx

shiftingScrIdx = find(shiftLabelsScr == 1);
noShiftingScrIdx = find(shiftLabelsScr == 0);

%after shuffling the data, recompute the average distance to the nearest
%five neighbors - same as above

% preallocate

nbShiftScr = cell(numel(shiftingScrIdx),1);
edShiftScr = cell(numel(shiftingScrIdx),1);

%compute the neighbors!

for n = 1:numel(shiftingScrIdx)
    [nbShiftScr{n},edShiftScr{n}] = knnsearch(proj_dataMat(shiftingScrIdx,:),proj_dataMat(shiftingScrIdx(n),:),'k',6);    
end

%do the same for no shifting!

nbNoShiftScr = cell(numel(noShiftingScrIdx),1);
edNoShiftScr = cell(numel(noShiftingScrIdx),1);

for n = 1:numel(noShiftingScrIdx)
    [nbNoShiftScr{n},edNoShiftScr{n}] = knnsearch(proj_dataMat(noShiftingScrIdx,:),proj_dataMat(noShiftingScrIdx(n),:),'k',6);
end

%recompute the mean distances after the scrambling of the assignments

meanEdShiftScr = zeros(numel(shiftingScrIdx),1); %preallocate
meanEdNoShiftScr = zeros(numel(noShiftingScrIdx),1);

%first the shift scramble
for n = 1:length(edShiftScr)
    meanEdShiftScr(n) = mean(edShiftScr{n}(2:6));
end

%no shift scramble

for n = 1:length(edNoShiftScr)
    meanEdNoShiftScr(n) = mean(edNoShiftScr{n}(2:6));
end


return

