function [r,log_llh_test,log_llh_test_model,log_llh_test_mean] = LLH_Comp(param,shift_A,spikes)

%% Description
% This code will take the LLH of the test data given the learned parameters
% from a previous A matrix


%% Initialize matrices and section the data for k-fold cross-validation

  
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
r = exp(shift_A * param'); n = spikes'; meanFR_test = nanmean(spikes); 
    
log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
log_llh_test = (-log_llh_test_model + log_llh_test_mean);
log_llh_test = log(2)*log_llh_test;
    
    
