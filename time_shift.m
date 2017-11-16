function [Shift_LLH_Max,p_shift,param,total_LLH_1] = time_shift(spiketrain,A,num_var,matComps,test_ind,train_ind,bestModelFits,param,variables)
%% convert A from cell array to matrix

A = cell2mat(A);

%create matrix form of matComps to enabling calling of only A matrix
%components of interest
matCompsMatForm = cell2mat(matComps);

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

%for all of the variables, saves the LLH values from the re-learning steps
%from each cross validation
total_ten_LLH_1 = zeros(numFolds,num_var);

%for all variables, saves the LLH values from the intial, unshifted set
%(and eventually, preceding re-learnign steps) from each cross validation
total_ten_LLH_0 = bestModelFits;

%saves the initial ten LLH values from the unshifted matrix (this will
%also be set to have the LLH values from preceding re-learning steps)
ten_LLH_0 = bestModelFits(:,end);

%this one is when we are testing all of the different shifts and computing
%the LLH
%this loop specifically places the different number of shifts into the
%first column so we can compare the LLH between different shifts
total_LLH_1 = zeros(numel(shift_vec),2,numel(shift_num));
for i = 1:size(total_LLH_1,3)
    total_LLH_1(:,1,i) = shift_vec;
end


%% Step 1: find the shift which maximizes the LLH of the fit

%this whole component will be iterated for each variable of interest within
%the A matrix
%preallocate A_shift
A_shift = A;

%preallocate for the p value from the signrank test
p_shift = cell(1,num_var);

% for all of variables
for k = 1:num_var 
    if variables(k) == 4  %since we don't want to shift the theta, if the variable being
        %shifted is theta, we will not shift it
        A_shift(:,matComps{k}) = A(:,matComps{k}); %don't shift it
        total_LLH_1(:,1,k) = zeros(length(shift_vec),1); %fill this in with zeros since we dont shift
        fprintf('\t\t- Not shifting - this variable is THETA\n');
    else
        %for all of the elements in the shift vector
        for m = 1:numel(shift_vec)
        %shifts that component of the matrix by the amount called from the
        %shift vector
        
        
            A_shift(:,matComps{k}) = circshift(A(:,matComps{k}),shift_vec(m));
            fprintf('\t\t- Shift %d\n', shift_vec(m));
        

            %compute the LLH of the fit, calling only the components of
            %interest with matCompsMatForm
            [~,log_llh_test,~,~] = LLH_Comp(param,A_shift(:,matCompsMatForm),spiketrain);

            %save all these data in the appropriate spot
            total_LLH_1(m,2,k) = log_llh_test;
        end
    end
    
    %now, find the maximum LLH - first, index where the maximum occurs
    [~,MaxLLH_Idx] = max(total_LLH_1(:,2,k));
    
    %returns the maximum
    Shift_LLH_Max(k) = total_LLH_1(MaxLLH_Idx,1,k);
    
    %reshift the A_shift to reflect the shift which maximized the fit
    
    %again, if variables == 4, dont do this
    if variables(k) == 4
        A_shift(:,matComps{k}) = A(:,matComps{k});
        fprintf('\t\t- Not re-shifting - variable is THETA');
    else
        A_shift(:,matComps{k}) = circshift(A(:,matComps{k}),Shift_LLH_Max(k));
    end
    %% Relearn the parameters, given the A_shift
    %again, we are only indexing the A matrix of interest
    [testFit,~,param1] =  fit_model_kfold(A_shift(:,matCompsMatForm),spiketrain,test_ind,train_ind);
    
    %finds the new LLH values from the k cross fold validation from the
    %second model fit
    ten_LLH_1 = testFit;
    total_ten_LLH_1(:,k) = ten_LLH_1;
    total_ten_LLH_0(:,k) = ten_LLH_0;
       
    %determine whether there is a significant difference in the LLH
    [p_shift{k}] = signrank(ten_LLH_1,ten_LLH_0);
    
    %if the LLH from the second fit is significantly larger than the first
    %fit
    if mean(ten_LLH_1) > mean(ten_LLH_0) && p_shift{k} < .05
        param = param1;
        A(:,matCompsMatForm) = A_shift(:,matCompsMatForm);
        ten_LLH_0 = ten_LLH_1;
        fprintf('The learned parameters SIGNIFICANTLY INCREASED the LLH of Fit.\n');
    else
        fprintf('The learned parameters did not significantly increase the LLH of Fit.\n');
      
    end 
   
end
    
return
    
     
    
    
    
    
    
    
    



