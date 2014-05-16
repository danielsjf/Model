%% FUNCTION SCEN  
function [scen_lim,scen_red] = SCEN(n_scen_gen,MU,eps_init,coef_fit_var,norm_var,forecast,samples_scenario,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,n_scen_select,actual_realization,ff_mode)   
    tic
    %% Scenario generation
    disp('    Scenario generation...')
    
    %n_scen_gen,MU,eps_init,coef_fit_var,norm_var,forecast,samples_scenario,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,n_scen_select,actual_realization,ff_mode

    % Generating histogram (without the function)
    forecast_pb = zeros(samples_scenario,1); % forecast per bin
    ll_int = [0:w_int:1-w_int]; % lower limits for the bins
    ul_int = [0+w_int:w_int:1]; % upper limits for the bins
    n_int = 1/w_int; % number of bins
    for b = 1:samples_scenario
        xx = find(forecast(b) >= ll_int & forecast(b) < ul_int);
        if xx >= n_int
            xx = n_int; 
        elseif xx <= 1
            xx = 1; 
        end
        forecast_pb(b) = xx;
   end

    % Exponential estimation of the correlation matrix
    x = mean(forecast_pb);
    eps = eps_init*sum(coef_fit_var.*[x^2 x 1])/norm_var;
    ETA = zeros(samples_scenario);
    for a = 1:samples_scenario
    for b = 1:samples_scenario
        ETA(a,b) = exp(-abs(a-b)/eps);
    end
    end
    
    % Generate random number
    RAND = mvnrnd(MU,ETA,n_scen_gen);

    % Translate to WPFE - scenarios
    inv_transform_matrix = zeros(n_scen_gen,samples_scenario);
    scen_error_pb = zeros(n_scen_gen,samples_scenario);
    scen_error = zeros(samples_scenario,n_scen_gen);

    for a = 1:n_scen_gen
       for b = 1:samples_scenario
          inv_transform_matrix(a,b) = normcdf(RAND(a,b));
          scen_error_pb(a,b) = find(cdf_stable(forecast_pb(b),:)>= inv_transform_matrix(a,b),1,'first');
          scen_error(b,a) = error_vec(scen_error_pb(a,b))+w_int/2;
       end   
    end  
    
    % add the actual realization and a zero-error scenario
    scen_error = [actual_realization zeros(samples_scenario,1) scen_error];
    n_scen_gen = n_scen_gen+2;    

    %% Probability Calculation
    disp('    Probability calculation...')
    
    % Initalisation 
    scen_prob_mat = zeros(samples_scenario*2,n_scen_gen);

    % Assign a probability to each element in a scenario
    for b = 1:n_scen_gen
       count = 0;
       for a = 1:samples_scenario
          % Probability of an error, given the forecast
          % Find the error bin the error belongs to
          error_pb = ceil((scen_error(a,b)+1)/w_int);
          % Look-up the probability
          scen_prob_mat(a+count,b) = pdf_stable(forecast_pb(a),error_pb);

          % Probability of an error, given the previous error
          % Find the error bin of the previous error and the error bin of the
          % current error
          if a == 1
              error_pb_e_prev = ceil((0+1)/w_ee);
          else
              error_pb_e_prev = ceil((scen_error(a-1,b)+1)/w_ee); 
          end      
          error_pb_e = ceil((scen_error(a,b)+1)/w_ee); 
          scen_prob_mat(a+count+1,b) = pdf_stable_ee(error_pb_e_prev,error_pb_e); 

          % position
          count = count+1;
       end   
    end

    % Probability of a scenario
    scen_prob = prod(scen_prob_mat);
    
    % Correct for impossible scenarios
    for b = 1:n_scen_gen
        for a = 1:samples_scenario
        if forecast(a)+scen_error(a,b) < 0 % negative production
        scen_prob(b) = 0;
        end
        if forecast(a)+scen_error(a,b) > 1 % production higher than installed capacity
            scen_prob(b) = 0;
        end
        end
    end
    
%     % Eliminate the scenarios with a probability equal to zero
%     del = find(scen_prob == 0);
%     if isempty(del) == 0
%     if del(1) == 1
%         del(1) = [];
%     end
%     scen_prob(del) = [];
%     scen_error(:,del) = [];
%     n_scen_gen = length(scen_prob);
%     end
    
%     % normalization
%     tot_prob = sum(scen_prob);
%     scen_prob = scen_prob./tot_prob;
    
   % equiprobability
    scen_prob = 1/n_scen_gen*ones(n_scen_gen,1)';
    
    % remove the actual realization and the zero-error scenario
    remove_error = scen_error(:,1:2);
    scen_error(:,1:2) = [];
    remove_prob = scen_prob(1:2);
    scen_prob(1:2) = [];
    n_scen_gen  = n_scen_gen - 2; 
    
 
%     % normalization
%     tot_prob = sum(scen_prob);
%     scen_prob = scen_prob./tot_prob;
   
    % Retain the scenarios with a probability above 10^-10
    cut_off_prob = 0; % probability limit to retain scenarios
    n_scen_gen_retain = 100; % nmber of generated scenarios that will be retained

    [useless,ind_scen] = find(scen_prob>cut_off_prob,n_scen_gen_retain,'first');
    scen_error = scen_error(:,ind_scen);
    
    % equiprobability
    n_scen_gen  = size(scen_error,2);
    scen_prob = 1/n_scen_gen*ones(n_scen_gen,1)';
    
%     scen_prob = scen_prob(ind_scen);%/sum(scen_prob(ind_scen));
    n_scen_gen = length(scen_prob);
            
    %% Scenario reduction 
    disp('    Scenario reduction: Fast Forward...')
    
    % Initialize the set of selected scenarios
    select_scen = zeros(n_scen_select,1);

    % initial cost matrix original scenario
    cost_mat = zeros(n_scen_gen,n_scen_gen);
    for jj=1:samples_scenario  
    cost_mat_t = abs(repmat(scen_error(jj,:),n_scen_gen,1)- repmat(scen_error(jj,:)',1,n_scen_gen));
    cost_mat = cost_mat+cost_mat_t;
    end
    cost_mat_orig = cost_mat;
    scen_prob_orig = scen_prob;

    if strcmp(ff_mode,'Advanced')==1
    % embedded variability scenarios
    scen_var = zeros(size(scen_error));
    for ii = 2:samples_scenario
       scen_var(ii,:) = abs(scen_error(ii,:)-scen_error(ii-1,:));
    end

    % cost matrix for variability
    cost_var = zeros(n_scen_gen,n_scen_gen); % cost function per time period
    for jj=1:samples_scenario  
    cost_mat_var= abs(repmat(scen_var(jj,:),n_scen_gen,1)- repmat(scen_var(jj,:)',1,n_scen_gen));
    cost_var = cost_var+cost_mat_var;
    end

    % added cost of variability
    cost_mat = cost_mat - cost_var; 
    % equiprobability
    scen_prob = ones(1,length(scen_prob_orig))./length(scen_prob_orig);
    end

    % Calculate the Kantorovich distance
    kantorovich=scen_prob*cost_mat; % kantorovich distance     

    % Select the first scenario that minimizes the Kantorovich distance
    [useless,select_scen(1)] = min(kantorovich); 

    for kk = 2:n_scen_select
        % Update the cost matrix: compare to the cost vector of the previously
        % selected scenario.
        selected = repmat(cost_mat(:,select_scen(kk-1)),1,n_scen_gen);    
        cost_mat_new = min(cost_mat,selected);
        cost_mat_new(select_scen(kk-1),:) = cost_mat(select_scen(kk-1),:);
        cost_mat = cost_mat_new;

        % the set of deleted scenarios
        del_scen = [1:1:n_scen_gen]'; 
        del_scen(select_scen(1:kk-1)) = [];

        % Calculate the Kantorovich distance
        kantorovich_calc = scen_prob(del_scen)*cost_mat(del_scen,del_scen); 
        kantorovich = 10^10*ones(n_scen_gen,1);
        kantorovich(del_scen) = kantorovich_calc;

        % Select the next scenario with minimum Kantorovich distance
        [useless,select_scen(kk)] = min(kantorovich); 
    end

    % Reduced set of scenarios:
    scen_error_red = scen_error(:,select_scen);
    scen_prob_red = scen_prob_orig(:,select_scen);

    % Not selected scenarios: 
    del_scen = [1:1:n_scen_gen]'; 
    del_scen(select_scen) = [];

    % Update the probabilities: transfer the probabilities from the non selected scenarios to the
    % selected ones.
    for a = 1:n_scen_gen-n_scen_select
    % Find the closest selected scenario
    [useless,ind] = min(cost_mat_orig(del_scen(a),select_scen),[],2);
    % Add the probabilites of the deleted scenarios to the selected scenarios
    % that are closest
    scen_prob_red(ind) = scen_prob_red(ind)+scen_prob_orig(del_scen(a));
    end
    scen_prob = scen_prob_orig;
    
    
    %% Add the zero-error-scenario to the reduced set of scenarios if needed
    scen_error_red = [remove_error(:,2) scen_error_red];
    scen_prob_red = [remove_prob(2) scen_prob_red];
    scen_prob_red = scen_prob_red/sum(scen_prob_red);
    
    %% Add the actual realization to the limited set of scenarios
    scen_error = [remove_error(:,1) scen_error];
    scen_prob = [remove_prob(1) scen_prob];
    scen_prob = scen_prob/sum(scen_prob);
    
    time = toc;
    disp(['    Done! Time elapsed: ',num2str(time),' s'])
  
    % Output
    scen_lim.error = scen_error;
    scen_lim.prob = scen_prob;
    
    scen_red.error = scen_error_red;
    scen_red.prob = scen_prob_red;
end