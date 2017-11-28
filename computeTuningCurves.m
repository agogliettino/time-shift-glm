% computes tuning curves for cells from time shift analysis

%% load appropriate data

load('cellData_11_27_17Run.mat'); %loads the data into the workspace - this contains all of the cells

%load the cells of interest
load('encodingCells.mat');

%load the variable grids
load('varGrids.mat');

%load the variables vectors
load('varVecs.mat');

%% recompute the design matrix

%since these variables were computed only on the first run (since we are
%concerned only with the dimensions of A, and are therefore treating it is
%a generic variable), we are going to use data from cell # 4 ie the first
%cell

A_TC = cell(4,1); %four because there are four variables
A_TC{1} = Allposgrid{4};
A_TC{2} = Allhdgrid{4};
A_TC{3} = Allspeedgrid{4};
A_TC{4} = Allthetagrid{4};

%% recompute the control points

%similar to above, control points (comprised of the vectors) are generic
%and do not differ between cells. so, we will use the fourth cell only

ctl_pts_all = cell(4,1); %four again because there are four variables

ctl_pts_all{1} = AllposVec{4};
ctl_pts_all{2} = AllhdVec{4};
ctl_pts_all{3} = AllspdVec{4};
ctl_pts_all{4} = AllthetaVec{4};

%% s for splines
s = 0.5;

%time in each bin
dt = 0.02;

%set plotfig to one to plot the curves
plotfig = 1;
%% plot tuning curves..
%call on the script using the above and then the variables_all from the
%glm

%must loop through all cells

% [tuning_curves] = plot_tuning(A,variables,parameters,ctl_pts_all,s,plotfig,dt)

    
[tuning_curves] = plot_tuning(A_TC,variables_all{4},AllParamVar{4},ctl_pts_all,s,1,dt);


