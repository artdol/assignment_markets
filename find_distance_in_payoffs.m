function [utils,closest_match,dist,closest_p]=find_distance_in_payoffs(x_point,mode,core_cons,A,c,ineq_mode)
% Finds the closest feasible point to a given payoff vector

[num_sellers, num_buyers] = size(A);

% utility as real-valued variables
u=sdpvar(1,num_buyers);
v=sdpvar(1,num_sellers);

% price vector as a vector of real-valued variables
p=sdpvar(1,num_sellers);

% all subsets of buyers
subsets_of_buyers = dec2bin(0:2^3-1) - '0';

% all permutations of goods
all_perms = perms([1 2 3]);

% start with a large distance as a candidate solution
min_dist = 10000;

for matched_buyers = subsets_of_buyers'
    % all permutations (assignments) of the goods to the subset of buyers
    for permutation_of_sellers = all_perms'
        
        % this generates every possible complete or incomplete matching
        % with slight overhead (some matchings may appear multiple times)
        mat = permutation_of_sellers;
        mat(matched_buyers == 0) = 0;
        
        % Add the feasibility constraints
        Cons=[];
        for i=1:num_sellers
            for j=1:num_buyers
                if mat(j)==i
                    if strcmp(ineq_mode,"weak") %"weak"
                        Cons=[Cons (u(j)+v(i)<=A(i,j)):char(strcat(string(i),string(j)))];
                    elseif strcmp(ineq_mode,"strong") %"strong"
                        Cons=[Cons (u(j)+v(i)==A(i,j)):char(strcat(string(i),string(j)))];
                    elseif strcmp(ineq_mode,"proj") %"projection"
                        Cons=[Cons (u(j)+v(i)<=A(i,j)):char(strcat(string(i),string(j)))];
                        r=sdpvar(1,1);
                        Cons=[Cons [v u]==x_point.*r, r>=0, r<=1];
                    end
                else
                     if core_cons % if true, add the core IC constraints
                        Cons=[Cons (u(j)+v(i)>=A(i,j))];
                     end
                end
            end
        end
        
        % unmatched players have zero payoff
        for j=1:num_buyers
            if mat(j)==0
                Cons=[Cons u(j)==0];
            end
        end
        Cons=[Cons v(setdiff(1:num_sellers,mat))==0];
        
        % prices are defined through costs and sellers' payoffs
        Cons=[Cons p==v+c];

        % optimization parameters. For cplex replace gurobi with cplex (or other
        % solvers)
        ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','bmibnb','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000);
        distcomp.feature( 'LocalUseMpiexec', false );

        % All variables are positive
        Cons=[Cons v>=0 u>=0 p>=0];

        if strcmp(mode,"euc")
        % Euclidian
        obj=norm(x_point-[v u]);

        elseif strcmp(mode,"chess")
        % Chessboard chebyshev
        obj=max(abs(x_point-[v u]));

        elseif strcmp(mode,"city")
        % Manhattan cityblock
        obj=sum(abs(x_point-[v u]));
        end

        if strcmp(ineq_mode,"proj")
            obj=-r;
        end

        % Find the closest feasible point
        solv=optimize([Cons],obj,ops);

        % Report if the program is not feasible
        if solv.problem~=0
            disp('Error: infeasible problem');
        end

        % check if the distance is less or equal the previous candidate
        if value(norm(x_point-[v u])) <= min_dist
            % if yes, update the candidate solution
            utils=value([v u]);
            dist=value(norm(x_point-[v u]));
            closest_match=mat;
            closest_p=value(p);
            min_dist = value(norm(x_point-[v u]));
        end
    end
end

end