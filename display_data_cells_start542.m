
%% load data

load('cellData_start542.mat'); %loads the data into the workspace

%load the cells of interest
load('encodingCells.mat');
%get rid of the crappy cells
encodingCells = encodingCells(encodingCells >= 542);

%find the indices of the significant p values for every shift

%preallocate for the significant p values and therefore shifts in the
%matrix

%% Find the significant shifts for each cell

SigAllVarShifts = AllVarShifts; %significant shifts

for i = encodingCells' %for all of the cells
   for j = numel(AllPShifts{i}) %for all the elements in the returned p vals
       if j > 0 %if there were any called shifts
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
        end
    end
end

%% Convert to matrix from cell array

posCells = cell2mat(posCells);
hdCells = cell2mat(hdCells);
spdCells = cell2mat(spdCells);

%% plot these data

figure(1);
bar(posCells(:,1),posCells(:,2),'FaceColor','k');
hold on
bar(hdCells(:,1),hdCells(:,2),'FaceColor','m');
hold on
bar(spdCells(:,1),spdCells(:,2),'FaceColor','c');
legend('Pos','HD','Spd');
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%% examine relationships in cells between the different variables and their respective shifts\

%first, group cells based on the different possible models

posMod = cell(numel(encodingCells),2); %1
hdMod = cell(numel(encodingCells),2); %2 
spdMod = cell(numel(encodingCells),2); %3
posHdMod = cell(numel(encodingCells),2); %1 2
hdSpdMod = cell(numel(encodingCells),2); %2 3
posSpdMod = cell(numel(encodingCells),2); %1 3
posHdSpdMod = cell(numel(encodingCells),2); %1 2 3

%create arrays for the multiple models 
ph = [1 2];
hs = [2 3];
ps = [1 3];
phs = [1 2 3];

for i = 1:numel(encodingCells)
    for j = 1:numel(variables_allFilt{i,2})
        if variables_allFilt{i,2}(j) == 4
            variables_allFilt{i,2}(j) = []; %gets rid of the theta model because we dont want it
            SigAllVarShiftsFilt{i,2}(j) = []; %gets rid of the theta shift, because we also don't want it
        end
    end
    if numel(variables_allFilt{i,2}) == 1 %if there is one variable in the model, do the next line
        if variables_allFilt{i,2} == 1
            posMod{i,1} = variables_allFilt{i,1};
            posMod{i,2} = SigAllVarShiftsFilt{i,2};
        elseif variables_allFilt{i,2} == 2
            hdMod{i,1} = variables_allFilt{i,1};
            hdMod{i,2} = SigAllVarShiftsFilt{i,2};
        elseif variables_allFilt{i,2} == 3
            spdMod{i,1} = variables_allFilt{i,1};
            spdMod{i,2} = SigAllVarShiftsFilt{i,2};
        end
    elseif numel(variables_allFilt{i,2}) == 2 %if there are two variables in the model, run these lines
        if variables_allFilt{i,2} == ph
            posHdMod{i,1} = variables_allFilt{i,1};
            posHdMod{i,2} = SigAllVarShiftsFilt{i,2};
        elseif variables_allFilt{i,2} == hs
            hdSpdMod{i,1} = variables_allFilt{i,1};
            hdSpdMod{i,2} = SigAllVarShiftsFilt{i,2};
        elseif variables_allFilt{i,2} == ps
            posSpdMod{i,1} = variables_allFilt{i,1};
            posSpdMod{i,2} = SigAllVarShiftsFilt{i,2};
        end
    elseif numel(variables_allFilt{i,2}) == 3 %if there are three variables in the model, run these
        if variables_allFilt{i,2} == phs 
            posHdSpdMod{i,1} = variables_allFilt{i,1};
            posHdSpdMod{i,2} = SigAllVarShiftsFilt{i,2};
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

%%display the data

%first the position
figure(2);
bar(posMod(:,1),posMod(:,2));
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');


%next, HD
figure(3);
bar(hdMod(:,1),hdMod(:,2));
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%next, spd
figure(4);
bar(spdMod(:,1),spdMod(:,2));
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');


%next, posHd
figure(5);
%pos first
bar(posHdMod(:,1),posHdMod(:,2),'FaceColor','m');
hold on
%next, HD
bar(posHdMod(:,1),posHdMod(:,3),'FaceColor','g');
legend('Pos','HD');
hold off
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%next, hdPsd
figure(6);
%hd firt
bar(hdSpdMod(:,1),hdSpdMod(:,2),'FaceColor','r');
hold on
%next, the Spd
bar(hdSpdMod(:,1),hdSpdMod(:,3),'FaceColor','c');
legend('HD','Spd');
hold off
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%next, PosSpd
figure(7);
%first the pos
bar(posSpdMod(:,1),posSpdMod(:,2),'FaceColor','k');
hold on
%next the spd
bar(posSpdMod(:,1),posSpdMod(:,3),'FaceColor','b');
legend('Pos','Spd');
hold off
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');

%last, the posHdSpd
figure(8);
%first the pos
bar(posHdSpdMod(:,1),posHdSpdMod(:,2),'FaceColor','g');
hold on
%next the Hd
bar(posHdSpdMod(:,1),posHdSpdMod(:,3),'FaceColor','r');
%next the spd
hold on
bar(posHdSpdMod(:,1),posHdSpdMod(:,4),'FaceColor','c');
legend('Pos','HD','Spd');
hold off
xlabel('Cell Number');
ylabel('Shift in A Mat (Max LLH)');








        
    