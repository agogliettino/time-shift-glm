function [Shift_LLH_Max,p,param,param1,total_LLH_1,ten_LLH_0,mean_LLH_1,mean_LLH_0,ten_LLH_1] = fit_shifted_dataPosHDSpd_OneFitGLM(spikes,A_mat,num_var,matComps,testFit,m,param)


%first, since there is no shift, set the initial shift = 0 - this depends
%on the initial input argument A_idx
shift_num = zeros(1,num_var);

%set Shift_LLH_Max to NaN to start with - the dimensions are also dependent
%on the input of A_idx
Shift_LLH_Max = nan(1,num_var);

%defines the number of shifts to be done in the A_mat matrix
shift_vec = -50:2:50;

%numfolds for the k fold cross validation
numFolds = 10;

total_ten_LLH_1 = zeros(numFolds,num_var);

total_ten_LLH_0 = zeros(numFolds,num_var);

%this one is when we are testing all of the different shifts and computing
%the LLH
%this loop specifically places the different number of shifts into the
%first column so we can compare the LLH between different shifts
total_LLH_1 = zeros(numel(shift_vec),2,numel(shift_num));
for i = 1:size(total_LLH_1,3)
    total_LLH_1(:,1,i) = shift_vec;
end


%gives us the ten LLH vals from the k fold cross validation
ten_LLH_0 = testFit(:,m);

%compute the mean initial LLH value
mean_LLH_0 = mean(ten_LLH_0);
%% Step 1: find the shift which maximizes the LLH of the fit

%this whole component will be iterated for each variable of interest within
%the A matrix
%preallocate A_shift
A_shift = A_mat;

%preallocate for the p value from the signrank test
p = cell(1,num_var);


% for all of variables
for k = 1:num_var
    %for all of the elements in the shift vector
    for m = 1:numel(shift_vec)
        %shifts that component of the matrix by the amount called from the
        %shift vector
        A_shift(:,matComps{k}) = circshift(A_mat(:,matComps{k}),shift_vec(m));
        fprintf('\t\t- Shift %d\n', shift_vec(m));
        
        %compute the LLH of the fit
        [~,log_llh_test,~,~] = LLH_Comp(param,A_shift,spikes);
        
        %save all these data in the appropriate spot
        total_LLH_1(m,2,k) = log_llh_test;
    end
    
    %now, find the maximum LLH - first, index where the maximum occurs
    [MaxLLH,MaxLLH_Idx] = max(total_LLH_1(:,2,k));
    
    %returns the maximum
    Shift_LLH_Max(k) = total_LLH_1(MaxLLH_Idx,1,k);
    
    %reshift the A_shift to reflect the shift which maximized the fit
    A_shift(:,matComps{k}) = circshift(A_mat(:,matComps{k}),Shift_LLH_Max(k));
    
    %% Relearn the parameters, given the A_shift
    [testFit,~,param1,~,~,~,~,~] = fit_model_kfold(A_shift,spikes,test_ind, train_ind);
    
    %finds the new LLH values from the k cross fold validation from the
    %second model fit
    ten_LLH_1 = testFit(:,3);
    
    total_ten_LLH_1(:,k) = ten_LLH_1;
    total_ten_LLH_0(:,k) = ten_LLH_0;
    
    %compute the mean
    mean_LLH_1 = mean(ten_LLH_1);
    
    %determine whether there is a significant difference in the LLH
    [p{k}] = signrank(ten_LLH_1,ten_LLH_0,'tail','right');
    
    %if the LLH from the second fit is significantly larger than the first
    %fit
    if mean_LLH_1 > mean_LLH_0 && p{k} < .05
        param = param1;
        A_mat = A_shift;
        mean_LLH_0 = mean_LLH_1;
        ten_LLH_0 = ten_LLH_1;
        fprintf('The learned parameters SIGNIFICANTLY INCREASED the LLH of Fit.\n');
    else
        fprintf('The learned parameters did not significantly increase the LLH of Fit.\n');
      
    end 
   
end
    

%% Plot the results

for i = 1:num_var
    figure(i)
    plot(shift_vec,total_LLH_1(:,2,i))
    xlabel('Shift in A matrix')
    ylabel('LLH')
end



end
    
     
    
    
    
    
    
    
    



