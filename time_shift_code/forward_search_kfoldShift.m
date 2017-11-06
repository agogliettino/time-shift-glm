function [allModelFits, bestModels, bestModelFits, parameters, pvals, final_pval,AllMaxLLHShifts] = forward_search_kfoldShift(A,spiketrain,train_ind,test_ind)

% do the forward search

num_var = length(A); % number of variables to search over

% do the forward-search method to identify variables encoded:

% var_vec = list of variables to add in
% variables = the variables currently in the model
% num_var = total number of variables to search over
% allModelFits = the model fit for every variable, on every iteration
% through the forward search procedure
% bestModelFits = the model fit for the best model (for every iteration)
% parametesr = parameters of the best model for each iteration
% pvals = the pvalue of the model comparison for each model

numFolds = length(test_ind);

% initialize values for while loop
var_vec = 1:num_var;
variables = [];
baseModel = -5*ones(numFolds,1);
allModelFits = [];
bestModelFits= [];
bestModels = [];
parameters = {};
pval = 0;
pvals = [];
loop_number = 0;
AllMaxLLHShifts = [];

while pval < 0.05 && numel(variables) < num_var
    
    testFit = nan(numFolds,num_var);
    trainFit = nan(numFolds,num_var);
    param_mean = cell(num_var,1);
    loop_number = loop_number + 1;

    for m = var_vec 

        
        fprintf('Fitting model  %d \n', m);
        
        % create matrix of variables in model currently
        X  = ones(length(spiketrain),1);
        temp_var = [variables m]; temp_var = sort(temp_var);
        for l = temp_var % this is always in order, for simplicity 
            X = [X A{l}];
        end
        
        [testFit(:,m),trainFit(:,m),param_mean{m}] = fit_model_kfold(X,spiketrain,test_ind,train_ind);
        
        %now, we have to shift the A matrix, to see if it increases the LLH
        %from the first fit, ie unshifted
        
        %set the input args according to which model it is
        
        %specifically sets the dimensions of the A matrix, depending on the
        %loop
        matComps = {size(A{m},2)};
        
        %call the function to shift the A mat
        [Shift_LLH_Max,~,~,~,total_LLH_1,~,~,~,~] = fit_shifted_dataPosHDSpd_OneFitGLM(spiketrain(test_ind,:),X(:,test_ind),num_var,matComps,testFit(:,m),m,param_mean{m});
        
        %determine if a different time bin (shift) maximized the fit
        %compared to unshifted
        
        %if there was a different shift, set the testFit from the previous
        %model fit to the new fit
        if Shift_LLH_Max ~=0
            testFit(:,m) = total_LLH_1;
        end
        
        % save the max Val shifts
        AllMaxLLHShifts(end+1) = Shift_LLH_Max;   
               
    end
       
        
        % save all of the model fits
        allModelFits = [allModelFits mean(testFit)'];
    
        % choose the best model
        [~,topModel_ind] = max(mean(testFit));
        topModel = testFit(:,topModel_ind);
        
        pval = signrank(topModel,baseModel,'tail','right');
        
        if pval < 0.05
            bestModelFits = [bestModelFits topModel ];
            bestModels = [bestModels topModel_ind ];
            parameters{end+1} = param_mean{topModel_ind};
        end
        
        pvals = [pvals pval];
        
        % find the variables in the model so far
        variables = [variables topModel_ind]; variables = sort(variables);

        % find the new variables to try adding into the model
        var_vec = setdiff(1:num_var,variables);

        baseModel = topModel;

end

% check that the final model is sig better than zero
final_pval = signrank(topModel,zeros(size(topModel)),'tail','right');

return




