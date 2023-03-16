function [blocking_cnt,blocking_cnt_traded]=check_PO_and_Blocks(ind_buy_pay,ind_sel_pay,present_s,A)
% Returns the blocking pairs for traded goods and overall

blocking_cnt=0;
blocking_cnt_traded=0;
for i=1:3
    for j=1:3
        if ind_buy_pay(j)+ind_sel_pay(i) <= A(i,j) - 1e-4
            blocking_cnt=blocking_cnt+1;
            if present_s(i)==1
                blocking_cnt_traded=blocking_cnt_traded+1;
            end
        end
    end
end

