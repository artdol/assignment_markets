%%% Parameters %%%

% Number of iterations
max_iter=50000;

% Matrix A
A=[[5 8 2];[7 9 6];[2 3 0]]*20;

% Game parameters
h=[[23 26 20];[22 24 21];[21 22 17]]*20;
c=[18,15,19]*20;

% bigM constant
bigM=10000;

% Maximum number of periods
T = 1500;

% Points for distance calculation
% leximin
lexi_p=[410 360 410];
% buyer-optimum 
buy_p=[420 400 380];
% seller-optimum
sell_p=[460   420   400];

% types of distances to calculate as an array
% possible modes are "euc" - Euclidean, "chess" - Chess distance, "city" -
% Manhattan cityblock distance. For example, to run all three use
% dist_modes={"euc","chess","city"};
% The first element will be  used for tables.
% The results in the paper use Euclidean distance.
dist_modes={"euc"};
num_dist_modes=length(dist_modes);

% core matching
core_m=get_nash_matching(A,[1 1 1],[1 1 1]);

% titles for types of agents to use as headings in the tables
agent_types = ["Uniform", "Highest margin", "Best reply"];

%%% Simulation %%%

% generate data
disp("Generating data (1/4)...")
[final_match,final_prices,u_vec,v_vec] = generate_simulated_data(max_iter,h,c,T);

% initialize solvers for closest point in the core to a given vector
% there is a solver for each distance regime
solv = cell(1,3);
for dist_ind=1:num_dist_modes
    dist_mode=dist_modes{dist_ind};
    solv{dist_ind}=init_Nash_distance_solver(dist_mode,A,h,c,core_m,[1 1 1],[1 1 1]);
end

% calculate distances to core and core refinements
disp("Calculating distances (core,2/4)...")
dist_to_core = zeros(3,max_iter,3);
dist_to_buy = zeros(3,max_iter,3);
dist_to_sell = zeros(3,max_iter,3);
dist_to_lex = zeros(3,max_iter,3);
for agent_type=1:3
    parfor iter=1:max_iter
        for dist_ind=1:num_dist_modes
            dist_mode=dist_modes{dist_ind};
            dist_to_core(agent_type,iter,dist_ind) = use_distance_solver(solv{dist_ind},double(final_prices(agent_type,:,iter)),dist_mode);
            dist_to_buy(agent_type,iter,dist_ind)=find_distance_to_point(double(final_prices(agent_type,:,iter)),buy_p,dist_mode);
            dist_to_sell(agent_type,iter,dist_ind)=find_distance_to_point(double(final_prices(agent_type,:,iter)),sell_p,dist_mode);
            dist_to_lex(agent_type,iter,dist_ind)=find_distance_to_point(double(final_prices(agent_type,:,iter)),lexi_p,dist_mode);
        end
    end
end


% calculate distances to Nash equilibria, L-core
disp("Calculating distances (Nash&Lcore,3/4)...")
dist_to_nash=zeros(3,max_iter,3);
dist_to_Lcore=zeros(3,max_iter,3);
for dist_ind=1:num_dist_modes
     dist_mode=dist_modes{dist_ind};
     
    % An array of solvers, one solver for each matching
    solv_nash = cell(4,4,4);
    solv_lcore = cell(4,4,4);
    
    for agent_type=1:3
        for iter=1:max_iter

            % print progress every 100 iterations
            if mod(iter,100)==0
                fprintf('%s(%d/3): %d of %d \n',dist_mode,agent_type,iter,max_iter);
            end

            % convert the matching to a cell array for use as an argument list
            matching=mat2cell(final_match(agent_type,:,iter)+1,1,ones(1,3));

            % initialize a new solver if none exists yet for this matching
            if isempty(solv_nash{matching{:}})
                
                present_b = [1 1 1]; % all buyers are present;
                present_s = ismember([1 2 3],final_match(agent_type,:,iter)); %active sellers = 1, inactive sellers = 0
            
                % find the Nash matching
                nash_m=get_nash_matching(A,present_b,present_s);
                
                % find the minimum excess, maximized in the L-core
                min_e=find_lcore(nash_m,present_s,A);

                % Initialize solvers for this set of sellers
                solv_nash{matching{:}}=init_Nash_distance_solver(dist_mode,A,h,c,nash_m,present_b,present_s);
                solv_lcore{matching{:}} = init_lcore_distance_solver(min_e,nash_m,present_s,dist_mode,A,h,c);
            end

            % use the solver to find the closest point to Nash eq.
            dist_to_nash(agent_type,iter,dist_ind) = use_distance_solver(solv_nash{matching{:}},squeeze(final_prices(agent_type,:,iter)),dist_mode);
            % find also the distance to L-core
            dist_to_Lcore(agent_type,iter,dist_ind)=find_distance_to_lcore(squeeze(final_prices(agent_type,:,iter)),solv_lcore{matching{:}},dist_mode);
        end
    end
 end
 
% Check the PO and blocking conditions
disp("Calculating the remaining theories and printing the results (4/4)...")

% Prepopulate results
blocking_cnt=zeros(3,max_iter,'uint8');
blocking_cnt_traded=zeros(3,max_iter,'uint8');
is_PO=zeros(3,max_iter,'uint8');
is_lottery_PO=zeros(3,max_iter,'uint8');

% Initialize Pareto Optimality solvers
solv = init_check_PO_solver(A);
solv_lott = init_check_lottery_PO_solver(A);

% run the solvers for every datapoint
for agent_type=1:3
    for iter=1:max_iter
        % print progress every 100 iterations
        if mod(iter,100)==0
            fprintf('%d/3: %d of %d \n',agent_type,iter,max_iter);
        end
        
        % Checking Pareto Optimality
        is_PO(agent_type,iter)=check_PO(u_vec(agent_type,:,iter),v_vec(agent_type,:,iter),solv);
        is_lottery_PO(agent_type,iter)=check_lottery_PO(u_vec(agent_type,:,iter),v_vec(agent_type,:,iter),solv_lott);
        
        % Counting blocking pairs
        [blocking_cnt(agent_type,iter),blocking_cnt_traded(agent_type,iter)]=check_Blocks(u_vec(agent_type,:,iter),v_vec(agent_type,:,iter),final_match(agent_type,:,iter),A);
    end
end

% Check if the points are Nash core for some weights
is_Ncore=zeros(3,max_iter);
zero_trades=zeros(3,max_iter);

% Initialize solvers that find Nash core weights, one solver for every full
% matching
full_matchings=perms([1 2 3]); 
solv=cell(1,length(full_matchings));

for agent_type=1:3
    for iter=1:max_iter
        matching=final_match(agent_type,:,iter);
        
        % Nash core can only rationalize full matchings and matchings where
        % the unmatched players have zero potential surplus
        if sum(matching>0)==3 || (all(matching == [1 2 0]) || all(matching == [2 1 0]))
            % In the latter case, match players for a zero surplus
            if (all(matching == [1 2 0]) || all(matching == [2 1 0]))
                m_corr=matching;
                m_corr(3)=3;
            else
                m_corr=matching;
            end
            
            % check if the solver exists, intialize it if not
            [~, ss_ind]=ismember(m_corr,full_matchings,'rows');
            if isempty(solv{ss_ind})
                solv{ss_ind}=init_Nash_core_solver(A,m_corr);
            end
            
            % run the solver
            is_Ncore(agent_type,iter)=check_Nash_core_weights(solv{ss_ind},double(u_vec(agent_type,:,iter)),double(v_vec(agent_type,:,iter)));
        end
        
        % count the number of zero-utility trades
        for b=1:3
            if u_vec(agent_type,b,iter)==0&&final_match(agent_type,b,iter)>0&&(b~=3||final_match(agent_type,b,iter)~=3)
                zero_trades(agent_type,iter)=zero_trades(agent_type,iter)+1;
            end
        end
        for s=1:3
            [~,matched_b] = ismember(s,final_match(agent_type,:,iter));
            if v_vec(agent_type,s,iter)==0&&ismember(s,final_match(agent_type,:,iter))&&(s~=3||matched_b~=3)
                zero_trades(agent_type,iter)=zero_trades(agent_type,iter)+1;
            end
        end
    end
end
        
% Find the Shapley value and Shapley values for subsets of markets
shval_array=zeros(6,7); % the Shapley value payoffs
shval_p_array=zeros(3,7,3); % the closest feasible price vector
for dist_ind=1:num_dist_modes 
    dist_mode=dist_modes{dist_ind};
    for i=1:7
        present_s=de2bi(i,3); % any given subset of sellers

        % find the Shapley value for this subset
        shval=get_shapley_value(A(present_s==1,1:3),sum(present_s)+3);

        % record the payoffs for the current subset of traded goods
        shval_array(4:6,i)=shval((sum(present_s)+1):(sum(present_s)+3)); %buyers
        shval_array(present_s==1,i)=shval(1:(sum(present_s))); %sellers

        % find the closest feasible price vector
        [~,~,~,shval_p_array(:,i,dist_ind)]=find_distance_in_payoffs(shval_array(:,i)',dist_mode,false,A,c,'strong');
    end
end

% Find distances to the Myerson-Shapley value
dist_to_myerson = zeros(3,max_iter,3);
for dist_ind=1:num_dist_modes
    dist_mode=dist_modes{dist_ind};
    for agent_type=1:3
        for iter=1:max_iter
            % Find the sold goods (present sellers)
            matching=final_match(agent_type,:,iter);
            present_s=zeros(1,3);
            present_s(matching(matching>0))=1;

            if sum(present_s)==0
                % if no goods are sold, the distance is zero
                dist_to_myerson(agent_type,iter)=0;
            else
                % otherwise, find the distance to the corresponding Shapley
                % value
                shval=shval_p_array(:,bi2de(present_s));
                dist_to_myerson(agent_type,iter,dist_ind)= find_distance_to_point(double(final_prices(agent_type,:,iter)),shval',dist_mode);
            end
        end
    end
end

% Count CE, NE, PO, etc. outcomes
count_ce=zeros(1,3);
count_ne=zeros(1,3);
count_PO=zeros(1,3);
count_lott_PO=zeros(1,3);
count_Lcore=zeros(1,3);
count_Ncore=zeros(1,3);
count_zero_trades=zeros(1,3);
eff_lvl=zeros(1,3);
sel_shr=zeros(1,3);
for agent_type=1:3
    count_ce(agent_type)=sum(blocking_cnt(agent_type,:)==0); 
    count_ne(agent_type)=sum(blocking_cnt_traded(agent_type,:)==0); 
    count_PO(agent_type)=sum(is_PO(agent_type,:)==1); 
    count_lott_PO(agent_type)=sum(is_lottery_PO(agent_type,:)==1); 
    eff_lvl(agent_type)=mean(sum(u_vec(agent_type,:,:),2)+sum(v_vec(agent_type,:,:),2));
    sel_shr(agent_type)=mean(squeeze(sum(v_vec(agent_type,:,:),2))./squeeze(sum(u_vec(agent_type,:,:),2)+sum(v_vec(agent_type,:,:),2)));
    count_Lcore(agent_type)=sum(dist_to_Lcore(agent_type,:,1)<=0.5); %uses euc. distance
    count_Ncore(agent_type)=sum(is_Ncore(agent_type,:)==1); 
    count_zero_trades(agent_type)=sum(zero_trades(agent_type,:)>0);
end

% print a table of matchings and their frequencies
% first find all possible matchings across agent types
matches=[];
for agent_type=1:3
    matches = [matches;unique(squeeze(final_match(agent_type,:,:))','rows')];
end
matches=unique(matches,'rows');

% count the matchings for each agent type
match_count=zeros(3,size(matches,1));
for agent_type=1:3
    for i=1:size(matches,1)
        match_count(agent_type,i)=sum(sum(matches(i,:)==squeeze(final_match(agent_type,:,:))',2)==3);
    end
end

% sort and display summary table (Table 3)
fprintf("Table 3 \n\n\n")
% heading
fprintf("%s & %s & %s & %s \\\\ \n", "Assignment (%% of observations)", agent_types)
fprintf("\\hline \n")
% matchings
[sorted_pairings_cnt, pairings_order] = sortrows(match_count','descend');
sorted_pairings=matches(pairings_order,:);
for i=1:length(sorted_pairings)
    matching=squeeze(sorted_pairings(i,:))';
    pairing=strcat(string(matching(1)),string(matching(2)),string(matching(3)));
    fprintf("[%s] & %0.2f & %0.2f & %0.2f \\\\ \n",pairing,sorted_pairings_cnt(i,1)/max_iter*100,sorted_pairings_cnt(i,2)/max_iter*100,sorted_pairings_cnt(i,3)/max_iter*100);
end
fprintf("\\hline \n")
% rationalization by theories
fprintf("Pareto Optimal & %0.2f & %0.2f & %0.2f \\\\ \n",count_PO/max_iter*100);
fprintf("PO, not dominated by lotteries & %0.2f & %0.2f & %0.2f \\\\ \n",count_lott_PO/max_iter*100);
fprintf("Core = CE = 0-blocked & %0.2f & %0.2f & %0.2f \\\\ \n",count_ce/max_iter*100);
fprintf("Competitive prices for traded goods = NE & %0.2f & %0.2f & %0.2f \\\\ \n",count_ne/max_iter*100);
fprintf("Least core for traded goods & %0.4f & %0.4f & %0.4f \\\\ \n",count_Lcore/max_iter*100);
fprintf("Pairwise Nash Core for some weights & %0.2f & %0.2f & %0.2f \\\\ \n",count_Ncore/max_iter*100);
fprintf("\\hline \n")
% blocking pairs
fprintf("1-blocked & %0.2f & %0.2f & %0.2f \\\\ \n",sum(blocking_cnt==1,2)/max_iter*100);
fprintf("2-blocked & %0.2f & %0.2f & %0.2f \\\\ \n",sum(blocking_cnt==2,2)/max_iter*100);
fprintf("3-blocked & %0.2f & %0.2f & %0.2f \\\\ \n",sum(blocking_cnt==3,2)/max_iter*100);
fprintf("\\ge 4-blocked & %0.2f & %0.2f & %0.2f \\\\ \n",sum(blocking_cnt>=4,2)/max_iter*100);
fprintf("\\hline \n")
% statistics
fprintf("Mean efficiency (%% of max total payoff) & %0.2f & %0.2f & %0.2f \\\\ \n",eff_lvl/320*100);
fprintf("Mean seller's share of surplus (%% of max total payoff) & %0.2f & %0.2f & %0.2f \\\\ \n",sel_shr*100);
fprintf("Outcomes with 0-utility trades & %0.2f & %0.2f & %0.2f \\\\ \n",count_zero_trades/max_iter*100);
fprintf("\\hline \n")

% display table of distances (Tables 11-13)
fprintf("Tables 11-13 \n\n\n")
% heading
fprintf("%s & %s & %s & %s \\\\ \n", "Assignment ", agent_types)
fprintf("\\hline \n")
% distances
is_core_assignment = sum([3 1 2]==squeeze(final_match(agent_type,:,:))',2)==3;
is_lex_assignment = sum([1 3 2]==squeeze(final_match(agent_type,:,:))',2)==3;

% by assignment
fprintf("Optimal & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_nash(:,is_core_assignment,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_nash(:,is_core_assignment,1),0,2));
fprintf("Leximin & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_nash(:,is_lex_assignment,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_nash(:,is_lex_assignment,1),0,2));
fprintf("All & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_nash(:,:,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_nash(:,:,1),0,2));

% other theories
fprintf("Shapley-Myerson value & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_myerson(:,:,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_myerson(:,:,1),0,2));
fprintf("Leximin & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_lex(:,:,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_lex(:,:,1),0,2));
fprintf("L-core & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_Lcore(:,:,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_Lcore(:,:,1),0,2));

% core refinements
fprintf("Seller-optimum & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_sell(:,is_core_assignment,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_sell(:,is_core_assignment,1),0,2));
fprintf("Buyer-optimum & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_buy(:,is_core_assignment,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_buy(:,is_core_assignment,1),0,2));
fprintf("L-core & %0.2f & %0.2f & %0.2f \\\\ \n",mean(dist_to_Lcore(:,is_core_assignment,1),2));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",std(dist_to_Lcore(:,is_core_assignment,1),0,2));

% Save the results in a matlab file
% Solver objects cannot be saved, so clear them first
clear solv solv_lcore solv_lott solv_nash;
save("powersimulations")

disp("Finished.")