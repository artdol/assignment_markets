% Matrix A
A=[[5 8 2];[7 9 6];[2 3 0]]*20;

% Game parameters
h=[[23 26 20];[22 24 21];[21 22 17]]*20;
c=[18,15,19]*20;

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

% first round for analysis (dropping the first 5 rounds)
round_begin=6; 

% load raw experimental data
T = load("raw_oTree_data",'T');
T=T.T;

% extract unique subject codes
subjects = unique(T.participant_code);
num_subjects = length(subjects);

% extract unique group codes
groups=unique(T.total_group);
num_groups=length(groups);

% total number of rounds
num_rounds = 15;

% first non-discarded round (discard the first 5)
first_round = 6;

% core matching
core_m=get_nash_matching(A,[1 1 1],[1 1 1]);

% treatment names for tables
treatment_names = {'DA','MP','PT','CI'};

% numer of observations
num_obs = zeros(1,4);

% parse payoff, treatment, prices, and matching information
disp("Parsing experimental data (1/4)...")
u_vec = zeros(num_groups,num_rounds,3);
v_vec = zeros(num_groups,num_rounds,3);
final_match = zeros(num_groups,num_rounds,3);
final_prices = zeros(num_groups,num_rounds,3);
group_treatment = cell(num_groups,1);
for gr_ind=1:num_groups
    group=groups(gr_ind);
    
    % player ids of the group
    players=unique(T.participant_code(strcmp(T.total_group,group)));
    
    % treatment code
    if T.subsession_treatment(strcmp(T.total_group,group))=="double_auction"
         group_treatment{gr_ind}='DA';
         num_obs(1) = num_obs(1) + num_rounds - first_round + 1;
    elseif T.subsession_treatment(strcmp(T.total_group,group))=="telephone_pit" 
         group_treatment{gr_ind}='PT';
         num_obs(3) = num_obs(3) + num_rounds - first_round + 1;
    elseif T.subsession_treatment(strcmp(T.total_group,group))=="complete_info" 
         group_treatment{gr_ind}='CI';
         num_obs(4) = num_obs(4) + num_rounds - first_round + 1;
    else
         group_treatment{gr_ind}='MP';
         num_obs(2) = num_obs(2) + num_rounds - first_round + 1;
    end
    
    % parse information for each round
    for round_ind=1:num_rounds
        
        if strcmp(group_treatment{gr_ind},'MP') || strcmp(group_treatment{gr_ind},'DA')
            % Auction treatments: parse final asks and bids in round
            asks=T.group_group_asks(strcmp(T.total_group,group)& T.subsession_round_number==round_ind);
            bids=T.group_group_bids(strcmp(T.total_group,group)& T.subsession_round_number==round_ind);
            % convert from strings to vectors
            asks = regexprep(asks,'[[],]','');
            asks = strsplit(asks{1});
            bids = regexprep(bids,'[[],]','');
            bids = strsplit(bids{1});
            bids=cellfun(@str2num, bids);
            asks=cellfun(@str2num, asks);
            % record the highest of the two
            final_prices(gr_ind,round_ind,:) = max(asks,bids);
        else
            % Bargaining treatments: parse offers in round
            for p_ind=1:length(players)
                p=players(p_ind);
                pid=T.player_id_in_group(strcmp(T.participant_code,p));
                
                if strcmp(T.player_player_type(strcmp(T.participant_code,p)),"seller")
                    pl_id=T.player_seller_id(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind)+1;
                    winning_bid = T.player_winning_bid(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);
                    
                    if isnan(winning_bid)
                        % if there is no winning offer for the good, take
                        % the highest of the standing offers by seller or
                        % any buyer
                        offers_from=T.player_player_offers(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);
                        offers_to=T.player_offers_to_player(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);
                        offers_from = regexprep(offers_from,'[[],]','');
                        offers_from = strsplit(offers_from{1});
                        offers_to = regexprep(offers_to,'[[],]','');
                        offers_to = strsplit(offers_to{1});
                        offers_to=cellfun(@str2num, offers_to);
                        offers_from=cellfun(@str2num, offers_from);
                        
                        % Since in all treatments going below the cost is
                        % not allowed, the price is at least as high as the
                        % cost to the seller c(pl_id)
                        final_prices(gr_ind,round_ind,pl_id) = max([c(pl_id) max(max(offers_from, offers_to))]);
                    else
                        final_prices(gr_ind,round_ind,pl_id) = winning_bid;
                    end
                end
            end
        end
        % parse payoffs and matches for each player
        for p_ind=1:length(players)
             p=players(p_ind);
             pid=T.player_id_in_group(strcmp(T.participant_code,p));
             
             
             if strcmp(T.player_player_type(strcmp(T.participant_code,p)),"buyer")
                 %buyer
                 pl_id=T.player_buyer_id(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind)+1;
                 u_vec(gr_ind,round_ind,pl_id)=T.player_payout(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);

                 % record the matched good
                 obj_id=T.player_object_traded(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);
                 if isnan(obj_id)
                     obj_id = 0;
                 end
                 final_match(gr_ind,round_ind,pl_id)=obj_id;
             else
                 %seller
                 pl_id=T.player_seller_id(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind)+1;
                 v_vec(gr_ind,round_ind,pl_id)=T.player_payout(strcmp(T.participant_code,p) & T.subsession_round_number==round_ind);
             end
        end
    end
end

% initialize solvers for closest point in the core to a given vector
% there is a solver for each distance regime
solv = cell(1,3);
for dist_ind=1:num_dist_modes
    dist_mode=dist_modes{dist_ind};
    solv{dist_ind}=init_Nash_distance_solver(dist_mode,A,h,c,core_m,[1 1 1],[1 1 1]);
end

% calculate distances to core and core refinements
disp("Calculating distances (core,2/4)...")
dist_to_core = zeros(num_groups,num_rounds,3);
dist_to_buy = zeros(num_groups,num_rounds,3);
dist_to_sell = zeros(num_groups,num_rounds,3);
dist_to_lex = zeros(num_groups,num_rounds,3);
for gr_ind=1:num_groups
    for round_ind=1:num_rounds
        for dist_ind=1:1
            dist_mode=dist_modes{dist_ind};
            dist_to_core(gr_ind,round_ind,dist_ind) = use_distance_solver(solv{dist_ind},squeeze(final_prices(gr_ind,round_ind,:))',dist_mode);
            dist_to_buy(gr_ind,round_ind,dist_ind)=find_distance_to_point(squeeze(final_prices(gr_ind,round_ind,:))',buy_p,dist_mode);
            dist_to_sell(gr_ind,round_ind,dist_ind)=find_distance_to_point(squeeze(final_prices(gr_ind,round_ind,:))',sell_p,dist_mode);
            dist_to_lex(gr_ind,round_ind,dist_ind)=find_distance_to_point(squeeze(final_prices(gr_ind,round_ind,:))',lexi_p,dist_mode);
        end
    end
end


% calculate distances to Nash equilibria, L-core
disp("Calculating distances (Nash&Lcore,3/4)...")
dist_to_nash=zeros(num_groups,num_rounds,3);
dist_to_Lcore=zeros(num_groups,num_rounds,3);
for dist_ind=1:num_dist_modes
     dist_mode=dist_modes{dist_ind};
     
    % An array of solvers, one solver for each matching
    solv_nash = cell(4,4,4);
    solv_lcore = cell(4,4,4);
    
    for gr_ind=1:num_groups
        for round_ind=1:num_rounds

            % convert the matching to a cell array for use as an argument list
            matching=mat2cell(squeeze(final_match(gr_ind,round_ind,:))'+1,1,ones(1,3));

            % initialize a new solver if none exists yet for this matching
            if isempty(solv_nash{matching{:}})
                
                present_b = [1 1 1]; % all buyers are present;
                present_s = ismember([1 2 3],squeeze(final_match(gr_ind,round_ind,:))'); %active sellers = 1, inactive sellers = 0
            
                % find the Nash matching
                nash_m=get_nash_matching(A,present_b,present_s);
                
                % find the minimum excess, maximized in the L-core
                min_e=find_lcore(nash_m,present_s,A);

                % Initialize solvers for this set of sellers
                solv_nash{matching{:}}=init_Nash_distance_solver(dist_mode,A,h,c,nash_m,present_b,present_s);
                solv_lcore{matching{:}} = init_lcore_distance_solver(min_e,nash_m,present_s,dist_mode,A,h,c);
            end

            % use the solver to find the closest point to Nash eq.
            dist_to_nash(gr_ind,round_ind,dist_ind) = use_distance_solver(solv_nash{matching{:}},squeeze(final_prices(gr_ind,round_ind,:))',dist_mode);
            % find also the distance to L-core
            dist_to_Lcore(gr_ind,round_ind,dist_ind)=find_distance_to_lcore(squeeze(final_prices(gr_ind,round_ind,:))',solv_lcore{matching{:}},dist_mode);
        end
    end
 end

% Check the PO and blocking conditions
disp("Calculating the remaining theories and printing the results (4/4)...")

% Prepopulate results
blocking_cnt=zeros(num_groups,num_rounds,'uint8');
blocking_cnt_traded=zeros(num_groups,num_rounds,'uint8');
is_PO=zeros(num_groups,num_rounds,'uint8');
is_lottery_PO=zeros(num_groups,num_rounds,'uint8');

% Initialize Pareto Optimality solvers
solv = init_check_PO_solver(A);
solv_lott = init_check_lottery_PO_solver(A);

% run the solvers for every datapoint
for gr_ind=1:num_groups
    for round_ind=1:num_rounds

        % Checking Pareto Optimality
        is_PO(gr_ind,round_ind)=check_PO(squeeze(u_vec(gr_ind,round_ind,:))',squeeze(v_vec(gr_ind,round_ind,:))',solv);
        is_lottery_PO(gr_ind,round_ind)=check_lottery_PO(squeeze(u_vec(gr_ind,round_ind,:))',squeeze(v_vec(gr_ind,round_ind,:))',solv_lott);

        % Counting blocking pairs
        [blocking_cnt(gr_ind,round_ind),blocking_cnt_traded(gr_ind,round_ind)]=check_Blocks(squeeze(u_vec(gr_ind,round_ind,:))',squeeze(v_vec(gr_ind,round_ind,:))',squeeze(final_match(gr_ind,round_ind,:))',A);
    end
end

% Check if the points are Nash core for some weights
is_Ncore=zeros(gr_ind,round_ind);
zero_trades=zeros(gr_ind,round_ind);

% Initialize solvers that find Nash core weights, one solver for every full
% matching
full_matchings=perms([1 2 3]); 
solv=cell(1,length(full_matchings));
for gr_ind=1:num_groups
    for round_ind=1:num_rounds
        matching=squeeze(final_match(gr_ind,round_ind,:))';
        
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
            is_Ncore(gr_ind,round_ind)=check_Nash_core_weights(solv{ss_ind},squeeze(u_vec(gr_ind,round_ind,:))',squeeze(v_vec(gr_ind,round_ind,:))');
        end
        
        % count the number of zero-utility trades
        for b=1:3
            if u_vec(gr_ind,round_ind,b)==0&&final_match(gr_ind,round_ind,b)>0&&(b~=3||final_match(gr_ind,round_ind,b)~=3)
                zero_trades(gr_ind,round_ind)=zero_trades(gr_ind,round_ind)+1;
            end
        end
        for s=1:3
            [~,matched_b] = ismember(s,squeeze(final_match(gr_ind,round_ind,:))');
            if v_vec(gr_ind,round_ind,s)==0&&ismember(s,squeeze(final_match(gr_ind,round_ind,:))')&&(s~=3||matched_b~=3)
                zero_trades(gr_ind,round_ind)=zero_trades(gr_ind,round_ind)+1;
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
dist_to_myerson = zeros(num_groups,num_rounds,3);
for dist_ind=1:num_dist_modes 
    dist_mode=dist_modes{dist_ind};
    for gr_ind=1:num_groups
        for round_ind=1:num_rounds
            % Find the sold goods (present sellers)
            matching=final_match(gr_ind,round_ind,:);
            present_s=zeros(1,3);
            present_s(matching(matching>0))=1;

            if sum(present_s)==0
                % if no goods are sold, the distance is zero
                dist_to_myerson(gr_ind,round_ind,dist_ind)=0;
            else
                % otherwise, find the distance to the corresponding Shapley
                % value
                shval=shval_p_array(:,bi2de(present_s));
                dist_to_myerson(gr_ind,round_ind,dist_ind)= find_distance_to_point(squeeze(final_prices(gr_ind,round_ind,:))',shval',dist_mode);
            end
        end
    end
end

% Count CE, NE, PO, etc. outcomes
count_ce=zeros(1,4);
count_ne=zeros(1,4);
count_PO=zeros(1,4);
count_lott_PO=zeros(1,4);
count_Lcore=zeros(1,4);
count_Ncore=zeros(1,4);
count_zero_trades=zeros(1,4);
eff_lvl=zeros(1,4);
sel_shr=zeros(1,4);
total_welfare_array=cell(1,4);
seller_welfare_array=cell(1,4);
count_blocked=zeros(4,4);
for treatment_ind=1:4
    for gr_ind=1:num_groups
            if strcmp(group_treatment{gr_ind},treatment_names{treatment_ind})
                count_ce(treatment_ind)=count_ce(treatment_ind)+sum(blocking_cnt(gr_ind,first_round:num_rounds)==0); 
                count_ne(treatment_ind)=count_ne(treatment_ind)+sum(blocking_cnt_traded(gr_ind,first_round:num_rounds)==0); 
                count_PO(treatment_ind)=count_PO(treatment_ind)+sum(is_PO(gr_ind,first_round:num_rounds)==1); 
                count_lott_PO(treatment_ind)=count_lott_PO(treatment_ind)+sum(is_lottery_PO(gr_ind,first_round:num_rounds)==1); 
                total_welfare_array{treatment_ind}=[total_welfare_array{treatment_ind} sum(u_vec(gr_ind,first_round:num_rounds,:),3)+sum(v_vec(gr_ind,first_round:num_rounds,:),3)];
                seller_welfare_array{treatment_ind}=[seller_welfare_array{treatment_ind} squeeze(sum(v_vec(gr_ind,first_round:num_rounds,:),3))./squeeze(sum(u_vec(gr_ind,first_round:num_rounds,:),3)+sum(v_vec(gr_ind,first_round:num_rounds,:),3))];
                count_Lcore(treatment_ind)=count_Lcore(treatment_ind)+sum(dist_to_Lcore(gr_ind,first_round:num_rounds,1)<=0.5); %uses euc. distance
                count_Ncore(treatment_ind)=count_Ncore(treatment_ind)+sum(is_Ncore(gr_ind,first_round:num_rounds)==1); 
                count_zero_trades(treatment_ind)=count_zero_trades(treatment_ind)+sum(zero_trades(gr_ind,first_round:num_rounds)>0);
                for num_blocks=1:3
                    count_blocked(treatment_ind,num_blocks)=count_blocked(treatment_ind,num_blocks)+sum(blocking_cnt(gr_ind,first_round:num_rounds)==num_blocks); 
                end
                count_blocked(treatment_ind,4)=count_blocked(treatment_ind,4)+sum(blocking_cnt(gr_ind,first_round:num_rounds)>=4); 

            end
    end
    eff_lvl(treatment_ind)=mean(total_welfare_array{treatment_ind});
    sel_shr(treatment_ind)=mean(seller_welfare_array{treatment_ind},'omitnan');
    
end

% print a table of matchings and their frequencies
% first find all possible matchings across agent types
matches=[];
for gr_ind=1:num_groups
    matches = [matches;unique(squeeze(final_match(gr_ind,first_round:num_rounds,:)),'rows')];
end
matches=unique(matches,'rows');

% count the matchings for each agent type
match_count=zeros(4,size(matches,1));
for treatment_ind=1:4
    for i=1:size(matches,1)
        for gr_ind=1:num_groups
            if strcmp(group_treatment{gr_ind},treatment_names{treatment_ind})
                match_count(treatment_ind,i)=match_count(treatment_ind,i)+sum(sum(matches(i,:)==squeeze(final_match(gr_ind,first_round:num_rounds,:)),2)==3);
            end
        end
    end
end

% sort and display summary table (Table 2)
fprintf("Table 2 \n\n\n")
% heading
fprintf("%s & %s & %s & %s & %s \\\\ \n", "Assignment (%% of observations)", treatment_names{:})
fprintf("\\hline \n")
% matchings
[sorted_pairings_cnt, pairings_order] = sortrows(match_count','descend');
sorted_pairings=matches(pairings_order,:);
for i=1:length(sorted_pairings)
    matching=squeeze(sorted_pairings(i,:))';
    pairing=strcat(string(matching(1)),string(matching(2)),string(matching(3)));
    fprintf("[%s] & %0.2f & %0.2f & %0.2f  & %0.2f \\\\ \n",pairing,sorted_pairings_cnt(i,1)/num_obs(1)*100,sorted_pairings_cnt(i,2)/num_obs(2)*100,sorted_pairings_cnt(i,3)/num_obs(3)*100,sorted_pairings_cnt(i,4)/num_obs(4)*100);
end
fprintf("\\hline \n")
% rationalization by theories
fprintf("Pareto Optimal & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_PO./num_obs*100);
fprintf("PO, not dominated by lotteries & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_lott_PO./num_obs*100);
fprintf("Core = CE = 0-blocked & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_ce./num_obs*100);
fprintf("Competitive prices for traded goods = NE & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_ne./num_obs*100);
fprintf("Least core for traded goods & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_Lcore./num_obs*100);
fprintf("Pairwise Nash Core for some weights & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_Ncore./num_obs*100);
fprintf("\\hline \n")
% blocking pairs
fprintf("1-blocked & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_blocked(:,1)'./num_obs*100);
fprintf("2-blocked & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_blocked(:,2)'./num_obs*100);
fprintf("3-blocked & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_blocked(:,3)'./num_obs*100);
fprintf("\\ge 4-blocked & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",count_blocked(:,4)'./num_obs*100);
fprintf("\\hline \n")
% statistics
fprintf("Mean efficiency (%% of max total payoff) & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",eff_lvl/320*100);
fprintf("Mean seller's share of surplus (%% of max total payoff) & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",sel_shr*100);
fprintf("Outcomes with 0-utility trades & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",(count_zero_trades)./num_obs*100);
fprintf("\\hline \n")


% organize distances and bargaining weights by treatments instead of groups
dist_to_nash_by_treatment = cell(4,3);
dist_to_myerson_by_treatment = cell(4,3);
dist_to_lex_by_treatment = cell(4,3);
dist_to_Lcore_by_treatment = cell(4,3);
dist_to_sell_by_treatment = cell(4,3);
dist_to_buy_by_treatment = cell(4,3);
dist_to_Lcore_given_optimal_assignment_by_treatment = cell(4,3);
dist_to_nash_given_optimal_assignment_by_treatment = cell(4,3);
dist_to_nash_given_leximin_assignment_by_treatment = cell(4,3);
for treatment_ind=1:4
    for dist_ind=1:num_dist_modes
        dist_mode=dist_modes{dist_ind};
        for i=1:size(matches,1)
            for gr_ind=1:num_groups
                if strcmp(group_treatment{gr_ind},treatment_names{treatment_ind})
                    for round_ind = first_round:num_rounds
                        dist_to_nash_by_treatment{treatment_ind,dist_ind} = [dist_to_nash_by_treatment{treatment_ind,dist_ind} dist_to_nash(gr_ind,round_ind,dist_ind)];
                        dist_to_myerson_by_treatment{treatment_ind,dist_ind} = [dist_to_myerson_by_treatment{treatment_ind,dist_ind} dist_to_myerson(gr_ind,round_ind,dist_ind)];
                        dist_to_lex_by_treatment{treatment_ind,dist_ind} = [dist_to_lex_by_treatment{treatment_ind,dist_ind} dist_to_lex(gr_ind,round_ind,dist_ind)];
                        dist_to_Lcore_by_treatment{treatment_ind,dist_ind} = [dist_to_Lcore_by_treatment{treatment_ind,dist_ind} dist_to_Lcore(gr_ind,round_ind,dist_ind)];
                        if all(squeeze(final_match(gr_ind,round_ind,:))' == [3 1 2])
                            dist_to_sell_by_treatment{treatment_ind,dist_ind} = [dist_to_sell_by_treatment{treatment_ind,dist_ind} dist_to_sell(gr_ind,round_ind,dist_ind)];
                            dist_to_buy_by_treatment{treatment_ind,dist_ind} = [dist_to_buy_by_treatment{treatment_ind,dist_ind} dist_to_buy(gr_ind,round_ind,dist_ind)];
                            dist_to_Lcore_given_optimal_assignment_by_treatment{treatment_ind,dist_ind} = [dist_to_Lcore_given_optimal_assignment_by_treatment{treatment_ind,dist_ind} dist_to_Lcore(gr_ind,round_ind,dist_ind)];
                            dist_to_nash_given_optimal_assignment_by_treatment{treatment_ind,dist_ind} = [dist_to_nash_given_optimal_assignment_by_treatment{treatment_ind,dist_ind} dist_to_nash(gr_ind,round_ind,dist_ind)];
                        end
                        if all(squeeze(final_match(gr_ind,round_ind,:))' == [ 1 3 2 ])

                            dist_to_nash_given_leximin_assignment_by_treatment{treatment_ind,dist_ind} = [dist_to_nash_given_leximin_assignment_by_treatment{treatment_ind,dist_ind} dist_to_nash(gr_ind,round_ind,dist_ind)];

                        end 
                    end
                end
            end
        end     
    end
end

% display table of distances (Tables H.10-H.12)
fprintf("Tables H.10-H.12 \n\n\n")
% heading
fprintf("%s & %s & %s & %s \\\\ \n", "Assignment ", treatment_names{:})
fprintf("\\hline \n")
% distances

% by assignment
fprintf("Optimal & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_nash_given_optimal_assignment_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_nash_given_optimal_assignment_by_treatment(:,1)'));
fprintf("Leximin & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_nash_given_leximin_assignment_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_nash_given_leximin_assignment_by_treatment(:,1)'));
fprintf("All & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_nash_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_nash_by_treatment(:,1)'));

% other theories
fprintf("Shapley-Myerson value & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_myerson_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_myerson_by_treatment(:,1)'));
fprintf("Leximin & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_lex_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_lex_by_treatment(:,1)'));
fprintf("L-core & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_Lcore_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_Lcore_by_treatment(:,1)'));

% core refinements
fprintf("Seller-optimum & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_sell_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_sell_by_treatment(:,1)'));
fprintf("Buyer-optimum & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_buy_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_buy_by_treatment(:,1)'));
fprintf("L-core & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n",cellfun(@(x) mean(x),dist_to_Lcore_given_optimal_assignment_by_treatment(:,1)'));
fprintf(" & (%0.2f) & (%0.2f) & (%0.2f) & (%0.2f) \\\\ \n",cellfun(@(x) std(x),dist_to_Lcore_by_treatment(:,1)'));

% Calculate the bargaining weights with minimal variance around the mean
% across periods
variance_of_weights=cell(1,4);
buyer_mean_weights=cell(1,4);
seller_mean_weights=cell(1,4);
buyer_weights=cell(1,4);
seller_weights=cell(1,4);

weights_dist_to_mean=cell(1,4);
weights_dist_to_sixth=cell(1,4);

for treatment_ind=1:4
    ind=1;
    flat_m=zeros(0,3);
    flat_u_vec=zeros(0,3);
    flat_v_vec=zeros(0,3);
    for gr_ind=1:num_groups
        if strcmp(group_treatment{gr_ind},treatment_names{treatment_ind})
            for round_ind=first_round:num_rounds
                matching=squeeze(final_match(gr_ind,round_ind,:))';
                if sum(matching>0)==3 || (all(matching == [1 2 0]) || all(matching == [2 1 0]))
                    if is_Ncore(gr_ind,round_ind)
                        if (all(matching == [1 2 0]) || all(matching == [2 1 0]))
                            m_corr=matching;
                            m_corr(3)=3;
                            
                            if u_vec(gr_ind,round_ind,3)>0  || v_vec(gr_ind,round_ind,3)>0
                                disp(gr_ind)
                                disp(round_ind)
                            end
                        else
                            m_corr=matching;
                        end
                        flat_m(ind,:)=m_corr;
                        flat_u_vec(ind,:)=u_vec(gr_ind,round_ind,:);
                        flat_v_vec(ind,:)=v_vec(gr_ind,round_ind,:);
                        ind=ind+1;
                        
                        
                    end
                end

            end
        end
    end
    num_points=ind-1;
    [variance_of_weights{treatment_ind},buyer_mean_weights{treatment_ind},seller_mean_weights{treatment_ind},buyer_weights{treatment_ind},seller_weights{treatment_ind}]=find_min_variance_bargaining_weights(num_points,flat_u_vec,flat_v_vec,A,flat_m);

    % record distances from the mean
    for p=1:num_points
        weights_dist_to_mean{treatment_ind}(p)=0;
        weights_dist_to_sixth{treatment_ind}(p)=0;
        for s=1:3
            weights_dist_to_mean{treatment_ind}(p)=weights_dist_to_mean{treatment_ind}(p)+(buyer_weights{treatment_ind}(p,b) - buyer_mean_weights{treatment_ind}(b))^2;
            weights_dist_to_sixth{treatment_ind}(p)=weights_dist_to_sixth{treatment_ind}(p)+(buyer_weights{treatment_ind}(p,b) - 1/6)^2;
        end
        for b=1:3
            weights_dist_to_mean{treatment_ind}(p)=weights_dist_to_mean{treatment_ind}(p)+(seller_weights{treatment_ind}(p,s) - seller_mean_weights{treatment_ind}(s))^2;
            weights_dist_to_sixth{treatment_ind}(p)=weights_dist_to_sixth{treatment_ind}(p)+(seller_weights{treatment_ind}(p,s) - 1/6)^2;
        end
    end
    % add the observations that are not Ncore rationalizable
    weights_dist_to_mean{treatment_ind}((num_points + 1) : (num_points + num_obs(treatment_ind) - count_Ncore(treatment_ind))) = 1;
    weights_dist_to_sixth{treatment_ind}((num_points + 1) : (num_points + num_obs(treatment_ind) - count_Ncore(treatment_ind))) = 1;

end

% calculate simple shares of payoffs
share_of_payoff_b = cell(4,3,3);
share_of_payoff_s = cell(4,3,3);
taking_turns = zeros(1,4);

share_is_50_50_flat = zeros(1,4);
share_is_50_50 = zeros(4,3,3);
share_is_50_50_210 = zeros(4,3,3);

matched_pairs = zeros(1,4);
matched_pairs_210 = zeros(1,4);
for treatment_ind=1:4
    for gr_ind=1:num_groups
        if strcmp(group_treatment{gr_ind},treatment_names{treatment_ind})
            for round_ind=first_round:num_rounds
                matching=squeeze(final_match(gr_ind,round_ind,:))';
                if (all(matching == [1 2 0]) || all(matching == [2 1 0]))
                    m_corr=matching;
                    m_corr(3)=3;
                else
                    m_corr=matching;
                end

                for b=1:3
                    s=m_corr(b);

                    if s>0
                        matched_pairs(treatment_ind) = matched_pairs(treatment_ind)+1;
                        if all(m_corr == [2 1 3]) 
                             matched_pairs_210(treatment_ind) = matched_pairs_210(treatment_ind)+1;
                        end
                        
                        if s~=3 || b~=3
                            % count episodes of subjects taking turns                    
                            if u_vec(gr_ind,round_ind,b) / A(s,b)== 1 || u_vec(gr_ind,round_ind,b) / A(s,b)== 0  
                                if u_vec(gr_ind,round_ind,b) / A(s,b)== 1- u_vec(gr_ind,round_ind-1,b) / A(s,b)
                                    taking_turns(treatment_ind) = taking_turns(treatment_ind)+1;
                                end
                            end

                           share_of_payoff_b{treatment_ind,b,s} = [share_of_payoff_b{treatment_ind,b,s} u_vec(gr_ind,round_ind,b) / A(s,b)];
                           share_of_payoff_s{treatment_ind,s,b} = [share_of_payoff_s{treatment_ind,s,b} v_vec(gr_ind,round_ind,s) / A(s,b)];

                           if u_vec(gr_ind,round_ind,b) / A(s,b) == .5
                                share_is_50_50_flat(treatment_ind) = share_is_50_50_flat(treatment_ind)+1;
                                share_is_50_50(treatment_ind,b,s) = share_is_50_50(treatment_ind,b,s) + 1;
                                if all(m_corr == [2 1 3]) 
                                  share_is_50_50_210(treatment_ind,b,s) = share_is_50_50_210(treatment_ind,b,s) + 1;
                                end
                           end
                        end
                    end
                end
            end
        end
    end
end

% make an array out of cell-array share_of_payoff
share_of_payoff_b_flat = cell(1,4);
share_of_payoff_s_flat = cell(1,4);

% confidence
share_of_payoff_b_flat_std = cell(1,4);
for treatment_ind=1:4
    for b=1:3
        for s=1:3
            share_of_payoff_b_flat{treatment_ind}(b,s)=mean(share_of_payoff_b{treatment_ind,b,s});
            share_of_payoff_b_flat_std{treatment_ind}(b,s)=std(share_of_payoff_b{treatment_ind,b,s});
            
            share_of_payoff_s_flat{treatment_ind}(s,b)=mean(share_of_payoff_s{treatment_ind,s,b});
        end
    end
    squeeze(share_is_50_50(treatment_ind,1,3)+share_is_50_50(treatment_ind,2,1))/sum(sum(share_is_50_50(treatment_ind,:,:)))
end


% plot the empricial CDFs of the distance of the empirical weights to the
% means by treatment
for treatment_ind=1:4
    cdfplot(weights_dist_to_mean{treatment_ind})
    hold on
end
hold off
legend(treatment_names)

% plot the empricial CDFs of the distance of the empirical weights to 
% 1/6 by treatment
for treatment_ind=1:4
    cdfplot(weights_dist_to_sixth{treatment_ind})
    hold on
end
hold off
legend(treatment_names)

% Save the results in a matlab file
% Solver objects cannot be saved, so clear them first
clear solv solv_lcore solv_lott solv_nash;
save("experimentaldata")

% Optionally: plot the spiderweb graph (Figure 3) of bargaining weights
% spider_plot function by Moses is available at https://de.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
for t=1:4
    fig = spider_plot([[buyer_weights{t},seller_weights{t}];[buyer_mean_weights{t},seller_mean_weights{t}]],...
    'Color', [repmat([0,0.5,0],length(buyer_weights{t}),1);[0.5,0,0]],...
    'AxesLabels', {'B1','B2','B3','S1','S2','S3'},...
    'AxesLimits', [0, 0, 0, 0, 0, 0; .7, .7, .7, .7, .7, .7], ...
    'LineTransparency',[repelem(0.3,length(buyer_weights{t})),1],'MarkerTransparency',[repelem(0.3,length(buyer_weights{t})),1]);
    
    % each spiderweb is saved as a separate svg image
    saveas(fig,strcat(treatment_names{t},".svg") ) 
end

disp("Finished.")
