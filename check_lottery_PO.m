function is_PO=check_lottery_PO_v2(ind_buy_pay,ind_sel_pay,solv)
% run the solver, if there is no solution then return false
is_PO=true;
[~,prob_code,~]=solv(ind_buy_pay,ind_sel_pay);
if prob_code==0
    is_PO=false;
    return
end

