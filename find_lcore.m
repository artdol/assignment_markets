function min_e=find_lcore(nash_m,present_s,A)
%finds the minimum excess in the L-core

% inactive sellers cannot create surplus
A(present_s==0,:)=0;

% the excess of players
e = sdpvar(1,6);

% utilities as real variables
u = sdpvar(1,3);
v = sdpvar(1,3);

% the minimum excess encoded as variable
z=sdpvar(1,1);

% Constraints defining players' excesses
Cons=[];
for i=1:3
    if present_s(i)==1
        for j=1:3
            if i~=nash_m(j)
                Cons=[Cons e(i) <= u(i) - max(A(i,j) - v(j),0)]; % excess of seller
                Cons=[Cons e(j+3) <= v(j) - max(A(i,j) - u(i),0)]; % excess of buyer
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
        Cons=[Cons e(j+3) <= v(j)+eps]; % excess of buyer
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

% z is defined to be the minimum excess of all players, which is
% implemented through inequalities
for k=1:3
    if nash_m(k)>0 
        Cons=[Cons z <= e(3+k)];
    end
end
for k=1:3
    if present_s(k)==1
        Cons=[Cons z <= e(k)];
    end
end

% All variables are positive
Cons=[Cons z>=0 u>=0 v>=0];

% Maximize the minimum excess
obj=-z;

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);
solv=optimize([Cons],obj,ops);

% return the minimum excess
min_e=-value(obj);