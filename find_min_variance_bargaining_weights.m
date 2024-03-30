function [vari,beta_mean_val,alpha_mean_val,beta_val,alpha_val]=find_min_variance_bargaining_weights(num_points,ind_buy_pay,ind_sel_pay,A,m)

% small positive epsilon
eps=0.0001;

% bargaining weights as real variables, one for each player and every
% datapoint
alpha=sdpvar(num_points,3,'full');
beta=sdpvar(num_points,3,'full');

% the mean weight for every player
alpha_mean = sdpvar(1,3);
beta_mean = sdpvar(1,3);

% binary matrix indicating whether the seller or the buyer precludes the match
x=binvar(num_points,3,3,'full');

% the big M constant
bigM=1000;

% the constraints define the Pairwise Nash core for every datapoint
Cons=[];
for p=1:num_points
    for b=1:3
        for s=1:3
            Cons=[Cons ind_buy_pay(p,b)*(alpha(p,s)+beta(p,b))+x(p,s,b)*bigM>=A(s,b)*beta(p,b)-eps];
            Cons=[Cons ind_sel_pay(p,s)*(alpha(p,s)+beta(p,b))+(1-x(p,s,b))*bigM>=A(s,b)*alpha(p,s)-eps];
        end
    end
    Cons=[Cons sum(alpha(p,:)+beta(p,:))<=1+eps];
    Cons=[Cons sum(alpha(p,:)+beta(p,:))>=1-eps];
end

for p=1:num_points
    for b=1:3
        Cons=[Cons ind_buy_pay(p,b)*(alpha(p,m(p,b))+beta(p,b))==beta(p,b)*A(m(p,b),b)];
        Cons=[Cons ind_sel_pay(p,m(p,b))*(alpha(p,m(p,b))+beta(p,b))==alpha(p,m(p,b))*A(m(p,b),b)];
    end
end

Cons=[Cons alpha_mean>=eps alpha_mean<=1 beta_mean>=eps beta_mean<=1-eps ];

% objective value is the variance of the bargaining weights
obj=0;
for p=1:num_points
     for s=1:3
         obj=obj+(alpha(p,s)-alpha_mean(s))^2 ;
     end
     for b=1:3
         obj=obj+(beta(p,b)-beta_mean(b))^2 ;
     end
end

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000,'gurobi.TuneTimeLimit',Inf);

% All weights are within (0,1) open interval
Cons=[Cons alpha>=eps beta>=eps alpha<=1-eps beta<=1-eps ];

% Run the minimization program
solv = optimize([Cons],obj,ops);

if solv.problem>0 % the problem should always be feasible
    disp('Error: Infeasible problem or an issue with the solver');
end
        
% output variables
vari=value(obj)/(num_points-1); % minimized variance
beta_mean_val=value(beta_mean); % mean bargaining weights for buyers
alpha_mean_val=value(alpha_mean);% mean bargaining weights for sellers
beta_val=value(beta); % bargaining weights for buyers
alpha_val=value(alpha); % bargaining weights for sellers
