% Generate random solutions to the weighted Nash bargaining problem. The
% problem is not convex, and the program finds a feasible solution that
% maximizes the NBS objective (See Appendix for discussion).

% number of data points to generate
max_iter = 5000;

% utility as real-valued variables
u=sdpvar(1,3);
v=sdpvar(1,3);

% total payoff
w=sum(u)+sum(v);

% game description
A=[[5 8 2];[7 9 6];[2 3 0]];
h=[[23 26 20];[22 24 21];[21 22 17]];
c=[18,15,19];

% set of full matchings
permlist=perms([1:3]);

is_core = zeros(1,max_iter);
matches = zeros(max_iter,3);

for iter=1:max_iter
    % print progress every 100 iterations
    if mod(iter,100)==0
        fprintf('%d of %d \n',iter,max_iter);
    end
    
    % generate a random vector of weights
    alphas=[0 rand(1,5) 1];
    alphas=diff(sort(alphas)); % this ensures that the weights are drawn uniformly
    v_w=alphas(1:3);
    u_w=alphas(4:6);
    
    % initialize the candiate solution to zero
    best_obj = 0;

    % solve, holding matching constant for every full matching
    for perm_ind=1:size(permlist,1)
        Cons=[];
        % every pair divides utility according to the weights
        for g=1:3
            Cons=[Cons u(g)+v(permlist(perm_ind,g))==A(permlist(perm_ind,g),g)];
            Cons=[Cons u_w(g)*v(permlist(perm_ind,g))==v_w(permlist(perm_ind,g))*u(g)];
        end
       
        % All variables positive
        for i=1:3
             Cons=[Cons v(i)>=0 u(i)>=0];
        end

        obj=0; % no optimization, only feasibility

        % optimization parameters. For cplex replace gurobi with cplex (or other
        % solvers)
        ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','bmibnb','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001);
        distcomp.feature( 'LocalUseMpiexec', false );
        
        % run the solver
        solv=optimize([Cons],obj,ops);

        if solv.problem>0 % the problem should always be feasible
            disp('Error: Infeasible problem or an issue with the solver');
        end

        % compare the solution to the current candidate in terms of the
        % Nash bargaining product
        if best_obj<value(prod(v)*prod(u))
            best_obj=value(prod(v)*prod(u));
            
            % record the matching
            matches(iter,:) = permlist(perm_ind,:);

            % check the core conditions
            is_core(iter)=true;
            for b=1:3
                for s=1:3
                    if A(s,b)>value(u(b)+v(s))
                        is_core(iter)=false;
                    end
                end
            end
        end
    end
end

% count the matchings
match_count=zeros(1,size(matches,1));

for i=1:size(permlist,1)
    match_count(i)=sum(sum(permlist(i,:)==matches,2)==3);
end

% sort and display the frequency table (Table 4)
% heading
fprintf("%s & %s \\\\ \n", "Assignment (\% of observations)", agent_types)
fprintf("\hline")
% matchings
[sorted_pairings_cnt, pairings_order] = sortrows(match_count','descend');
sorted_pairings=permlist(pairings_order,:);
for i=1:length(sorted_pairings)
    matching=squeeze(sorted_pairings(i,:))';
    pairing=strcat(string(matching(1)),string(matching(2)),string(matching(3)));
    fprintf("[%s] & %0.2f%% \\\\ \n",pairing,sorted_pairings_cnt(i)/max_iter*100);
end
fprintf("\hline")

% optionally produce a scatterplot of solutions
% scatter3(u_iter(:,1),u_iter(:,2),u_iter(:,3))