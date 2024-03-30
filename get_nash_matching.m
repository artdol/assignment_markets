function Nash_m = get_nash_matching(A,present_b,present_s)

X = binvar(size(A,1),size(A,2),'full'); % matching integer variables

% matching constraints
Cons=[sum(X,1)<=1 sum(X,2)<=1];

% missing players
for i=1:3
    if present_b(i)==0
        Cons=[Cons sum(X(:,i))==0];
    end
    if present_s(i)==0
        Cons=[Cons sum(X(i,:))==0];
    end
end

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);
distcomp.feature( 'LocalUseMpiexec', false );

% objective value = total surplus
obj = -sum(sum(X.*A));

% solve the optimization problem
solv=optimize([Cons],obj,ops);

if solv.problem > 0 
    disp("Could not solve the problem")
end

% return the matching
Nash_m = [0 0 0];
for i=1:3
    if isempty(find(value(X(:,i)), 1)) || present_s(find(value(X(:,i)), 1))==0
        Nash_m(i) = 0;
    else
        Nash_m(i) = find(value(X(:,i)));
    end
end
end