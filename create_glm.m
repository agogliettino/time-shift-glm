function [LLH,param,variables, best_shift, shift_pval] = create_glm(posfile,spikefile,filt_eeg,sampleRate,boxSize)

% load spike file
load(spikefile)
dt = 0.02;

%% compute the input matrix

%%%%%%%%% POSITION %%%%%%%%%

fprintf('Making position spline\n');

% load and rescale position file
[posx,posy,posx2,posy2,post] = rescalePos(posfile,boxSize,filt_eeg,sampleRate);
posx_c = nanmean([posx posx2],2); posy_c = nanmean([posy posy2],2); % compute average position

bin_p = 10; s = 0.5;
posVec = linspace(0,boxSize,bin_p); posVec(1) = -0.01;
[posgrid,~] = spline_2d(posx_c,posy_c,posVec,s);
A{1} = posgrid;

%%%%%%%%% HEAD DIRECTION %%%%%%%%%
fprintf('Making head direction spline\n');

% compute the head direction
direction = atan2(posy2-posy,posx2-posx)+pi/2;
direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers


bin_h = 15;
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

spdVec = [0:5:50 60]; spdVec(1) = -0.1;
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

bin_t = 15; s = 0.5;
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
testFit = bestModelFits(:,end);
param = parameters{end};

%%%%%%%% DO THE TIME SHIFT PROCEDURE %%%%%%%
A_k = {};
var_sort = sort(variables);
for j = 1:numel(var_sort)
    A_k{end+1} = A{var_sort(j)};
end

shift_var = find(var_sort <= 3);
[best_shift1,shift_pval1,param0,LLH0] = find_timeshift(A_k,shift_var,spiketrain,param,testFit,test_ind,train_ind);
best_shift = zeros(1,3); best_shift(var_sort(shift_var)) = best_shift1;
shift_pval = zeros(1,3); shift_pval(var_sort(shift_var)) = shift_pval1;

param = param0;
LLH = LLH0;

return