function solv = init_check_PO_solver(A)

% the Big M constant
bigM=1000;

[num_sellers, num_buyers] = size(A);

% utility as real-valued variables
u=sdpvar(1,num_sellers,'full');
v=sdpvar(1,num_buyers,'full');

% matching as a binary matrix
matrix_m=binvar(3,3,'full');

% Parameters, specifying the utility vector that is checked for PO
ind_sel_pay=sdpvar(1,num_sellers,'full');
ind_buy_pay=sdpvar(1,num_buyers,'full');

% Is improvement possible?
% All players receive at least the same utility
Cons=[];
for j=1:3
    Cons=[Cons v(j)>=ind_buy_pay(j)-1e-4];
end
for i=1:3
    Cons=[Cons u(i)>=ind_sel_pay(i)-1e-4];
end
for j=1:3
    Cons=[Cons v(j)<=1e-4 + sum(matrix_m(:,j))*bigM];
end
for i=1:3
    Cons=[Cons u(i)<=1e-4 + sum(matrix_m(i,:))*bigM];
end
for i=1:num_sellers
    for j=1:num_buyers
        Cons=[Cons u(i)+v(j)<=A(i,j)+1e-4+(1-matrix_m(i,j))*bigM];
    end
end

% Some player strictly improves
Cons=[Cons sum(u)+sum(v)>=sum(ind_buy_pay)+sum(ind_sel_pay)+0.01];

obj=0; % no optimization, only feasibility

% optimization parameters. For cplex replace gurobi with cplex (or other
% solvers)
ops = sdpsettings('verbose',0,'savesolveroutput',1,'solver','gurobi','gurobi.TimeLimit',60*60*24,'gurobi.MIPGap',0.000001,'bmibnb.maxiter',1000);

% Return the optimizer object
solv=optimizer([Cons],obj,ops,{matrix_m,ind_buy_pay,ind_sel_pay},{u v});

end