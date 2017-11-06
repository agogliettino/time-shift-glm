function [posx,posy,posx2,posy2,post] = rescalePos(posfile,boxSize,filt_eeg,sampleRate)

% load data
load(posfile);

%% if there isn't a second LED recording, reject all of this
if numel(posx2) < 1 || numel(posy2) < 1
    posx = NaN; posy = NaN; posx2 = NaN; posy2 = NaN;
    return
end

%% take out NaN's and replace them with neighboring values
% interpolate to replace NaN's
% this only works if the first and last values aren't nan's
if isnan(posx(1))
    posx(1) = posx(find(isfinite(posx),1,'first'));
end
if isnan(posy(1))
    posy(1) = posy(find(isfinite(posy),1,'first'));
end
if isnan(posx2(1))
    posx2(1) = posx2(find(isfinite(posx2),1,'first'));
end
if isnan(posy2(1))
    posy2(1) = posy2(find(isfinite(posy2),1,'first'));
end

if isnan(posx(end))
    posx(end) = posx(find(isfinite(posx),1,'last'));
end
if isnan(posy(end))
    posy(end) = posy(find(isfinite(posy),1,'last'));
end
if isnan(posx2(end))
    posx2(end) = posx2(find(isfinite(posx2),1,'last'));
end
if isnan(posy2(end))
    posy2(end) = posy2(find(isfinite(posy2),1,'last'));
end

posx(isnan(posx)) = interp1(post(~isnan(posx)),posx(~isnan(posx)),post(isnan(posx)),'linear');
posy(isnan(posy)) = interp1(post(~isnan(posy)),posy(~isnan(posy)),post(isnan(posy)),'linear');
posx2(isnan(posx2)) = interp1(post(~isnan(posx2)),posx2(~isnan(posx2)),post(isnan(posx2)),'linear');
posy2(isnan(posy2)) = interp1(post(~isnan(posy2)),posy2(~isnan(posy2)),post(isnan(posy2)),'linear');

%% make sure that the length of the vectors are correct given the eeg recording
maxTime = numel(filt_eeg)/sampleRate*50;
posx = posx(1:min(maxTime,numel(post))); posx2 = posx2(1:min(maxTime,numel(post)));
posy = posy(1:min(maxTime,numel(post))); posy2 = posy2(1:min(maxTime,numel(post)));
post = post(1:min(maxTime,numel(post)));

%% rescale the x and y position, so position goes from 0 to boxsize (typically 50-100 cm)
% assumes shape is a box
% rescale x-position
maxval_x = max(max([posx posx2])); minval_x = min(min([posx posx2]));
posx = boxSize * (posx-minval_x) / (maxval_x-minval_x); 
posx2 = boxSize * (posx2-minval_x) / (maxval_x-minval_x); 

% rescale y-position
maxval_y = max(max([posy posy2])); minval_y = min(min([posy posy2]));
posy = boxSize * (posy-minval_y) / (maxval_y-minval_y); 
posy2 = boxSize * (posy2-minval_y) / (maxval_y-minval_y); 


return