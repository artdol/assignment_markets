function solv=init_Nash_core_solver(A,m)
% epsilon close to zero
eps=0.0001;

% bargaining weights as optimization variables
alpha=sdpvar(1,3);
beta=sdpvar(1,3);

% Parameters, specifying the utility vector that is checked for Nash core
ind_buy_pay=sdpvar(1,3);
ind_sel_pay=sdpvar(1,3);

% binary matrix indicating whether the seller or the buyer precludes the match
x=binvar(3,3,'full');

% big M constant
bigM=100000;

Cons=[];
% for i=1:3
%     Cons=[Cons gamma(m(i))+1/N/2>=alpha(m(i))+beta(i)];
%     Cons=[Cons gamma(m(i))-1/N/2<=alpha(m(i))+beta(i)];
% end

for b=1:3
    for s=1:3
        Cons=[Cons ind_buy_pay(b)*(alpha(s)+beta(b))+x(s,b)*bigM>=A(s,b)*beta(b)-eps];
        Cons=[Cons ind_sel_pay(s)*(alpha(s)+beta(b))+(1-x(s,b))*bigM>=A(s,b)*alpha(s)-eps];
    end
end

for b=1:3
    Cons=[Cons ind_buy_pay(b)*(alpha(m(b))+beta(b))==beta(b)*A(m(b),b)];
    Cons=[Cons ind_sel_pay(m(b))*(alpha(m(b))+beta(b))==alpha(m(b))*A(m(b),b)];
end

% weights are between epsilon and 1
Cons=[Cons alpha>=eps*100 beta>=eps*100 alpha<=1-eps*100 beta<=1-eps*100];

obj=0; % No optimization, only feasibility

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);

% return the optimizer object
solv=optimizer([Cons],obj,ops,{ind_buy_pay ind_sel_pay},{alpha beta x});

