% a script to do a PCA on the data to identify clusters of neurons within
% each variable

%% load the data

%loads the cells
load('encodingCells.mat'); 

%load the data on which the glm was run
load('do_glm_output_112817.mat');

%% create n x m matrix

%for each variable, we will create an n x m matrix, where n is the number
%of cells and m is the number of parameters. for each cell, we will place
%the model-derived firing rate for each time bin

%% Find the significant shifts for each cell

SigAllVarShifts = AllVarShifts; %significant shifts

for i = encodingCells' %for all of the cells
   for j = numel(AllPShifts{i}) %for all the elements in the returned p vals
       if j > 0 %if there were any called shifts i e having more than zero elements
            if AllPShifts{i}{j} > .05 %if the p value is not significant
                SigAllVarShifts{i}(j) = 0; %set the shift to zero
                
            end
       end
   end
end

%% remove empty cells

SigAllVarShiftsFilt = cell(numel(encodingCells),2); %preallocates
for k = 1:numel(encodingCells)
    SigAllVarShiftsFilt{k,2} = SigAllVarShifts{encodingCells(k)}; %in the filtered array, this will only contain non-empty vals
    SigAllVarShiftsFilt{k,1} = encodingCells(k); %places the cell name there
end

%% Count up the shifts across all cells

%here, we will plot the distribution of shifts, across all cells, for each
%of the variables: Pos, HD, and Spd

%first, filter out the empty variables

variables_allFilt = cell(numel(encodingCells),2); %preallocate for filtered variables

for i = 1:numel(encodingCells)
    variables_allFilt{i,2} = variables_all{encodingCells(i)}; %this will contain only non-empty vals
    variables_allFilt{i,1} = encodingCells(i); %puts the cell name there
end

%% organize cells into three different categories: those with Pos, HD, and Spd variables, and therefore shift data

posCells = cell(numel(encodingCells),2); %preallocate
hdCells = cell(numel(encodingCells),2);
spdCells = cell(numel(encodingCells),2);
thetaCells = cell(numel(encodingCells),2);

%place the cell name and appropriate shift into the array

for f = 1:numel(encodingCells) %for all the number of cells
    for l = 1:numel(variables_allFilt{f,2}) %for all the elements in each variable cell
        if variables_allFilt{f,2}(l) == 1 %if it is Pos
            posCells{f,1} = SigAllVarShiftsFilt{f,1}; %put the cell name in the array
            posCells{f,2} = SigAllVarShiftsFilt{f,2}(l); %put the relevent shift in the array
        elseif variables_allFilt{f,2}(l) == 2 %if it is HD
            hdCells{f,1} = SigAllVarShiftsFilt{f,1}; %put the cell name into the array
            hdCells{f,2} = SigAllVarShiftsFilt{f,2}(l); %place relevant shifts in array
        elseif variables_allFilt{f,2}(l) == 3 %if it is Spd
            spdCells{f,1} = SigAllVarShiftsFilt{f,1}; %putt the cell name into array
            spdCells{f,2} = SigAllVarShiftsFilt{f,2}(l); %put the shift in array
        elseif variables_allFilt{f,2}(l) == 4 %if its theta
            thetaCells{f,1} = SigAllVarShiftsFilt{f,1};
            thetaCells{f,2} = SigAllVarShiftsFilt{f,2}(l);
        end
    end
end

%% Convert to matrix from cell array
posCells = cell2mat(posCells);
hdCells = cell2mat(hdCells);
spdCells = cell2mat(spdCells);
thetaCells = cell2mat(thetaCells);

boxSizes = xlsread('Wildtype_Mice_Database_v5.xlsx');
boxSizeCell = boxSizes(:,13);

cellboxSize100 = find(boxSizeCell == 100); %cells with box size of 100

variables_all100 = {};
for i = 1:length(cellboxSize100)
    variables_all100{end+1} = variables_all{i};
end

%find size of the tuning curve for Pos only
lengthIdxPos = []; 
for i = 1:length(variables_all100)
    if numel(variables_all100{i}) == 1 && variables_all100{i} == 1 %if it is a pos
        lengthIdxPos(end+1) = length(tuning_curves_all{i});
    end
end





%find position only tuning curves

%need to first find where the pos bins end

posOnly = [];

for i = encodingCells'
    if numel(variables_all{i}) == 1 && variables_all{i} == 1
        posOnly(end+1) = i;
    end
end




%% build pos matrix

%preallocate using largest num of cols
pMatPos = nan(length(posCells(:,1)),length(tuning_curves_all{lengthIdxTot(1)}));

%add the pos

%check the size
posSize = [];

for i = posOnly
    posSize(end+1) = size(tuning_curves_all{i},1);
end
    
        

for i = 1:numel(posCells(:,1))
    pMatPos(i,1:length(tuning_curves_all{posCells(i,1)})) = tuning_curves_all{posCells(i,1)};
end







