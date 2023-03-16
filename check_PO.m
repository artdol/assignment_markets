function is_PO=check_PO(ind_buy_pay,ind_sel_pay,solv)

% The set of all full matchings
all_m=perms([1 2 3]);
num_matchings=size(all_m,1);

% go over all matchings to find a possible improvement
for m=1:num_matchings
    cur_m=all_m(m,:);
    
    % convert a matching to the binary matrix form
    matrix_m=zeros(3);
    for j=1:3
        if cur_m(j)>0
            matrix_m(cur_m(j),j)=1;
        end
    end
    
    % run the solver, if there is no solution then return false
    is_PO=true;
    [~,prob_code,~]=solv(matrix_m,ind_buy_pay,ind_sel_pay);
    if prob_code==0
        is_PO=false;
        return
    end
end
