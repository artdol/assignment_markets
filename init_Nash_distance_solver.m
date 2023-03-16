function solv=init_core_Nash_distance_solver(dist_mode,A,h,c,nash_m,present_b,present_s)
% Initializes a solver for closest points in the core/Nash equilibrium set
% to any given vector.
% present_b, present_s, vectors of three, =1 if a buyer/seller is present
% in the problem, =0 if absent
% nash_m = the optimal matching to calculate the distance to
% dist_mode = distance regime

% utilities and prices as variables for optimization
u=sdpvar(1,3);
v=sdpvar(1,3);
p=sdpvar(1,3);

% the parameter vector, a point to which the distance is calculated
x_point=sdpvar(1,3);


% constraints defining utilities through prices
Cons=[];
for j=1:3
    for i=1:3
        if present_b(j)==1 && present_s(i)==1
            if nash_m(j)==i
                Cons=[Cons (v(j)==(h(i,j)-p(i)))];
                Cons=[Cons (u(i)==-c(i)+p(i))];
            end
        end
    end
end

% core constraints
for i=1:3
    for j=1:3
        if present_b(j)==1 && present_s(i)==1
            Cons=[Cons (u(i)+v(j)>=A(i,j)):char(strcat(string(i),string(j)))];
        end
    end
end

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000);
distcomp.feature( 'LocalUseMpiexec', false );

% all variables are positive
 for i=1:3
     Cons=[Cons v(i)>=0 u(i)>=0 p(i)>=0 ];
 end

% types of distances
if strcmp(dist_mode,"euc")
% Euclidian
obj=norm(x_point-p);

elseif strcmp(dist_mode,"chess")
% Chessboard chebyshev
obj=max(abs(x_point-p));

elseif strcmp(dist_mode,"city")
% Manhattan cityblock
obj=sum(abs(x_point-p));
end

% return the optimizer object
solv=optimizer([Cons],obj,ops,x_point,{u v p});

end