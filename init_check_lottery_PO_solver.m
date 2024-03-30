function solv = init_check_lottery_PO_solver(A)
[num_sellers, num_buyers] = size(A);

% The set of all full matchings
all_m=perms([1 2 3]);
num_matchings=size(all_m,1);

% utility as real-valued variables
u=sdpvar(num_matchings,num_sellers,'full');
v=sdpvar(num_matchings,num_buyers,'full');
% expected utility from lotteries
U=sdpvar(1,num_sellers,'full');
V=sdpvar(1,num_buyers,'full');

% lottery over matchings as a real variable
x=sdpvar(1,num_matchings,'full');

% Constraints define the utility as a function of matching 
Cons=[];
for m=1:num_matchings
    cur_m=all_m(m,:);

    for j=1:3
        if cur_m(j)>0
            Cons=[Cons u(m,cur_m(j))+v(m,j)==x(m)*A(cur_m(j),j)];
        else
            Cons=[Cons v(m,j)<=1e-4];
        end
    end
        
    Cons=[Cons u(m,setdiff(1:3,cur_m))==0];
end

Cons=[Cons sum(x)==1 x>=0]; % full matchings
Cons=[Cons u>=-1e-4 v>=-1e-4 U>=-1e-4 V>=-1e-4]; % positive utility
Cons=[Cons U==sum(u,1) V==sum(v,1)]; % expected utility

obj=-sum(U)-sum(V); % maximize the total payoff

% Parameters, specifying the utility vector that is checked for PO
ind_sel_pay=sdpvar(1,num_sellers,'full');
ind_buy_pay=sdpvar(1,num_buyers,'full');

% Is improvement over supplied utilities possible?
% All players receive at least the same utility
for j=1:3
    Cons=[Cons V(j)>=ind_buy_pay(j)-1e-4];
end
for i=1:3
    Cons=[Cons U(i)>=ind_sel_pay(i)-1e-4];
end

% Some player strictly improves
Cons=[Cons sum(U)+sum(V)>=sum(ind_buy_pay)+sum(ind_sel_pay)+0.01];

obj=0; % No optimization, only feasibility

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);

% Return the optimizer object
solv=optimizer([Cons],obj,ops,{ind_buy_pay,ind_sel_pay},{u v U V x});

end