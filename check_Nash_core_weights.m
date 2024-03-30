function [t,weights]=check_Nash_core_weights(solv,ind_buy_pay,ind_sel_pay)
% use the provided optimizer object
[weights,prob]=solv(ind_buy_pay,ind_sel_pay);

% return 1 if the solution exists, 0 if the problem is infeasible
t=prob==0;