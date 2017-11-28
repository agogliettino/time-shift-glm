function [param,tuning_curves,variables,VarShifts,P_ShiftVar,paramVar,VarTotal_LLH_1,ForwardPval,ForwardFinal_Pval,varAllModelFits,A,posgrid,hdgrid,speedgrid,thetagrid,posVec,hdVec,spdVec,thetaVec] = create_glm(posfile,spikefile,filt_eeg,sampleRate,boxSize)

% load spike file
load(spikefile)
dt = 0.02;

%% compute the input matrix

%%%%%%%%% POSITION %%%%%%%%%

fprintf('Making position spline\n');

% load and rescale position file
[posx,posy,posx2,posy2,post] = rescalePos(posfile,boxSize,filt_eeg,sampleRate);
posx_c = nanmean([posx posx2],2); posy_c = nanmean([posy posy2],2); % compute average position

bin_p = 8; s = 0.5;
posVec = linspace(0,boxSize,bin_p); posVec(1) = -0.01;
[posgrid,cpts_all] = spline_2d(posx_c,posy_c,posVec,s);
A{1} = posgrid;

%%%%%%%%% HEAD DIRECTION %%%%%%%%%
fprintf('Making head direction spline\n');

% compute the head direction
direction = atan2(posy2-posy,posx2-posx)+pi/2;
direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers


bin_h = 10;
hdVec = linspace(0,2*pi,bin_h+1); hdVec = hdVec(1:end-1);
s = 0.5;
[hdgrid] = spline_1d_circ(direction,hdVec,s);
A{2} = hdgrid;

%%%%%%%%% SPEED %%%%%%%%%
fprintf('Making speed spline\n');
% compute the speed
speed = nan(size(post));
for i = 1:numel(post)
    if i == 1
        speed(i) = sqrt((posx_c(i) - posx_c(i+1)).^2 + (posx_c(i) - posy_c(i+1)).^2)/(post(i+1) - post(i));
    elseif i == numel(post)
        speed(i) = sqrt((posx_c(i-1) - posx_c(i)).^2 + (posy_c(i-1) - posy_c(i)).^2)/(post(i) - post(i-1));
    else
        speed(i) = sqrt((posx_c(i-1) - posx_c(i+1)).^2 + (posy_c(i-1) - posy_c(i+1)).^2)/(post(i+1) - post(i-1));
    end
end
speed = conv(speed,gausswin(5)/sum(gausswin(5)),'same');

% set everything over 60 to 60
speed(speed > 60) = 60;

spdVec = [0:10:60]; spdVec(1) = -0.1;
s = 0.5;
[speedgrid,~] = spline_1d(speed,spdVec,s);

A{3} = speedgrid;

%%%%%%%%% THETA %%%%%%%%%
hilb_eeg = hilbert(filt_eeg); % compute hilbert transform
phase = atan2(imag(hilb_eeg),real(hilb_eeg)); %inverse tangent (-pi to pi)
ind = phase <0; phase(ind) = phase(ind)+2*pi; % from 0 to 2*pi

phase_ind = round(post*sampleRate);
phase_ind(phase_ind + 1>numel(filt_eeg)) = [];
theta_phase = phase(phase_ind+1); %gives phase of lfp at every time point

bin_t = 12; s = 0.5;
thetaVec = linspace(0,2*pi,bin_t+1); thetaVec = thetaVec(1:end-1);
thetagrid = spline_1d_circ(theta_phase,thetaVec,s);

A{4} = thetagrid;

%%%%%%%% SPIKE TRAIN %%%%%%%%%%%%%%%
% find # of spikes in every 20 ms bin
spiketrain = hist(cellTS,post);


%%%%%%% COMPUTE TEST AND TRAIN INDICES %%%%%
% calculate train and testing indices
% define the test and train indices
numFolds = 10; T = numel(post); numPts = round(1/dt); % 1 second
numSeg = ceil(T/(numFolds*numPts));
oneSeg = ones(numPts,1)*(1:numFolds);
new_ind = repmat(oneSeg(:),numSeg,1);
new_ind = new_ind(1:T);

test_ind = cell(numFolds,1);
train_ind = cell(numFolds,1);
for k = 1:numFolds
    test_ind{k} = find(new_ind == k);
    train_ind{k} = setdiff(1:T,test_ind{k});
end


%%%%%%%% FORWARD SEARCH PROCEDURE %%%%%%%%%

[allModelFits, bestModels, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,train_ind,test_ind);

%
variables = bestModels;
%sort them to be in order
testFit = bestModelFits(:,end);
param = parameters{end};

%%%%%%%% DO THE TIME SHIFT PROCEDURE %%%%%%%
%depending on the variables in the model, we have to put a column vector of
%ones before the first column in the A matrix
%sort variables because the A matrix is in the 1-4 variable order
variables = sort(variables);

%% A is changed here
orig_A = A; 
if variables(1) == 1
    A{1} = [ones(length(spiketrain),1) A{1}];
elseif variables(1) == 2
    A{2} = [ones(length(spiketrain),1) A{2}];
elseif variables(1) == 3
    A{3} = [ones(length(spiketrain),1) A{3}];
elseif variables(1) == 4
    A{4} = [ones(length(spiketrain),1) A{4}];
end

%compute the A indices - because the A matrix is in cell array form,
%we can compute them separately
A1Idx = 1:size(A{1},2);
A2Idx = A1Idx(end)+1:A1Idx(end)+size(A{2},2);
A3Idx = A2Idx(end)+1:A2Idx(end)+size(A{3},2);
A4Idx = A3Idx(end)+1:A3Idx(end)+size(A{4},2);

%concatenate them all together
All_AIdx = {A1Idx A2Idx A3Idx A4Idx};

%figure out what matComps should be

matComps = cell(1,numel(variables));
for i = 1:numel(variables)
    matComps{i} = All_AIdx{variables(i)};
end

%number of variables in the model
num_var = numel(variables);

%% Compute best shift
%only if there were variables returned

if numel(num_var) > 0
    [Shift_LLH_Max,p_shift,param,total_LLH_1] = time_shift(spiketrain,A,num_var,matComps,test_ind,train_ind,bestModelFits,param,variables);
end

VarShifts = Shift_LLH_Max;
P_ShiftVar = p_shift;
paramVar = param;
VarTotal_LLH_1 = total_LLH_1;
ForwardPval = pvals;
ForwardFinal_Pval = final_pval;
varAllModelFits = allModelFits;

ctl_pts_all{1} = posVec;
ctl_pts_all{2} = hdVec;
ctl_pts_all{3} = spdVec;
ctl_pts_all{4} = thetaVec;
plotfig = 0; % this should be set to 0
[tuning_curves] = plot_tuning(orig_A,variables,param,ctl_pts_all,s,plotfig,dt);


return