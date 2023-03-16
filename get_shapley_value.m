function phi=get_shapley_value(A,num_players)
% Find the Shapley value

[num_sellers,num_buyers] = size(A);

% utility as real-valued variables
u=sdpvar(1,num_sellers);
v=sdpvar(1,num_buyers);

% A matrix as a paramter
A_mat=sdpvar(num_sellers,num_buyers,'full');

% Incentive cosntraints
Cons=[];
for i=1:num_sellers
    for j=1:num_buyers
        Cons=[Cons u(i)+v(j)>=A_mat(i,j)];
    end
end

% Utility is positive
Cons=[Cons u>=0 v>=0];

% Minimize total payoff, this is the equivalent formulation to the core LP
obj=sum(u)+sum(v);

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000);

% Create the optimizer object
solv=optimizer([Cons],obj,ops,A_mat,[v,u]);

% calculate the phi value for every player
phi=zeros(1,num_players);
for i=1:num_players
    % the set of other players
    others=setdiff(1:num_players,i);
    % go over every coalition size
    for coal_size=1:(num_players-1)
        % all coalitions of this size without player i
        coal_wo_i=nchoosek(others,coal_size);
        for coal_i=1:size(coal_wo_i,1)
            % any given coalition without player i
            coal=coal_wo_i(coal_i,:);
            buyers=coal(coal>num_sellers)-num_sellers;
            sellers=coal(coal<=num_sellers);
            
            % The A matrix with the players from the coalition
            A_coal=zeros(num_sellers,num_buyers);
            A_coal(sellers,buyers)=A(sellers,buyers);
            
            % Maximum payoff of the coalition without i
            Mwoi=sum(solv(A_coal));
            
            % same coalition with player i
            coal=[coal_wo_i(coal_i,:) i];
            buyers=coal(coal>num_sellers)-num_sellers;
            sellers=coal(coal<=num_sellers);
            
            % The A matrix with the players from the coalition
            A_coal=zeros(num_sellers,num_buyers);
            A_coal(sellers,buyers)=A(sellers,buyers);

            % Maximum payoff of the coalition with i   
            Mwi=sum(solv(A_coal));
            
            % The Shapley value calculation
            phi(i)=phi(i)+(Mwi- Mwoi)*factorial(coal_size)*factorial(num_players-coal_size-1)/factorial(num_players);
        end
        
    end
end
