function [tuning_curves] = plot_tuning(A,variables,parameters,ctl_pts_all,s,plotfig,dt)

%% Description
% Given the variables, A (and thetagrid and shgrid), and the parameters,
% this will return the tuning curves for the cell
% if plotfig = 1, this will also plot the tuning curves

% NOTE: I just use A to compute the correct indexes

variables = sort(variables);
b0 = parameters(1);
param = parameters(2:end);

% position
total_ind = 0;
scale = zeros(5,1);
if ismember(1,variables)
    param_ind = size(A{1},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [pos_y,pos_x1,pos_x2] = spline_2d_plot(param1,ctl_pts_all{1},s);
    scale(1) = mean(mean(pos_y));
end

% head direction
if ismember(2,variables)
    param_ind = size(A{2},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [hd_y,hd_x] = spline_1d_circ_plot(param1,ctl_pts_all{2},s);
    scale(2) = mean(hd_y);
end

% speed
if ismember(3,variables)
    param_ind = size(A{3},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [spd_y,spd_x] = spline_1d_plot(param1,ctl_pts_all{3},s);
    scale(3) = mean(spd_y);
end

% theta phase
if ismember(4,variables)
    param_ind = size(A{4},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [theta_y,theta_x] = spline_1d_circ_plot(param1,ctl_pts_all{4},s);
    scale(4) = mean(theta_y);
end

tuning_curves = [];

if plotfig
    figure()
    
    subplot(1,5,1)
    if ismember(1,variables)
        scale_factor_ind = setdiff(variables,1); scale_factor = scale(scale_factor_ind);
        imagesc(exp(pos_y)*exp(b0)*prod(exp(scale_factor))/dt);
        temp = exp(pos_y)*exp(b0)*prod(exp(scale_factor));
        tuning_curves = [tuning_curves; temp(:)];
        colorbar; axis off
    end
    
    subplot(1,5,2)
    if ismember(2,variables)
        scale_factor_ind = setdiff(variables,2); scale_factor = scale(scale_factor_ind);
        plot(hd_x,exp(hd_y)*exp(b0)*prod(exp(scale_factor))/dt,'k','linewidth',2);
        box off
        xlabel('HD angle')
        ylabel('spikes/s')
        axis tight
        tuning_curves = [tuning_curves; exp(hd_y')*exp(b0)*prod(exp(scale_factor))];
    end
    
    subplot(1,5,3)
    if ismember(3,variables)
        scale_factor_ind = setdiff(variables,3); scale_factor = scale(scale_factor_ind);
        plot(spd_x,exp(spd_y)*exp(b0)*prod(exp(scale_factor))/dt,'k','linewidth',2);
        box off
        xlabel('speed (cm/s)')
        ylabel('spikes/s')
        axis tight
        tuning_curves = [tuning_curves; exp(spd_y')*exp(b0)*prod(exp(scale_factor))];
    end
    
    subplot(1,5,4)
    if ismember(4,variables)
        scale_factor_ind = setdiff(variables,4); scale_factor = scale(scale_factor_ind);
        plot(theta_x,exp(theta_y)*exp(b0)*prod(exp(scale_factor))/dt,'k','linewidth',2);
        box off
        xlabel('theta angle')
        ylabel('spikes/s')
        axis tight
        tuning_curves = [tuning_curves; exp(theta_y')*exp(b0)*prod(exp(scale_factor))];
    end

    
end

keyboard

return