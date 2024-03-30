function solv=init_lcore_distance_solver(min_e,nash_m,present_s,dist_mode,A,h,c)
% Creates a solver object for finding the closest L-core point to a given
% vector
% min_e is the precalculated minimum excess value, maximized in the L-core
% any point in the L-core should have at least this objective value

% Small constant fir strict inequalities
eps=0.0001;

% utilities as real variables
u=sdpvar(1,3);
v=sdpvar(1,3);

% prices as real variables
p=sdpvar(1,3);
% the excess of players
e=sdpvar(1,6); 

% the parameter vector, a point to which the distance is calculated
x_point=sdpvar(1,3);

% inactive sellers cannot create surplus
A(present_s==0,:)=0;
   
% Constraints defining players' excesses
Cons=[];
for i=1:3
    if present_s(i)==1
        for j=1:3
            if i~=nash_m(j)
                Cons=[Cons e(i) <= u(i) - max(A(i,j) - v(j),0)]; % seller excess
                Cons=[Cons e(j+3) <= v(j) - max(A(i,j) - u(i),0)]; % buyer excess
            else
                if nash_m(j)~=0
                    Cons=[Cons u(i) + v(j) == A(i,j)]; % core constraint
                end
            end
        end
        Cons=[Cons e(i) <= u(i)+eps];
    end
end
for j=1:3
    if nash_m(j)>0
        Cons=[Cons e(j+3) <= v(j)+eps];
    end
end

for k=1:3
    if nash_m(k)>0 
      Cons=[Cons e(k+3)>=min_e-eps];
    end
end
for k=1:3
    if present_s(k)==1
     Cons=[Cons e(k)>=min_e-eps];
    end
end

% unmatched players have zero utility
for j=1:3
    if nash_m(j)==0
        Cons=[Cons v(j)==0];
    end
end
for i=1:3
    if all(nash_m~=i)
        Cons=[Cons u(i)==0];
    end
end

% constraints defining utilities through prices  
for j=1:3
    for i=1:3
        if nash_m(j)==i
            Cons=[Cons (v(j)==(h(i,j)-p(i)))];
            Cons=[Cons (u(i)==-c(i)+p(i))];
        end
    end
end

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);
distcomp.feature( 'LocalUseMpiexec', false );

% All variables are positive
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
solv=optimizer([Cons],obj,ops,x_point,{v u p e});

