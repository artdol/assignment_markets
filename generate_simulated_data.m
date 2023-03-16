function [final_match,final_prices,u_vec,v_vec] = generate_simulated_data(max_iter,h,c,T)
% Generates data points according to "uniform", "highest margin", "best
% reply" strategies.

% Three agent simulations to run
agent_types = ["uniform", "highest margin", "best reply"];

% Prepopulate output matrices, int16 to save memory
final_match=zeros(3,3,max_iter,'uint16');
final_prices=zeros(3,3,max_iter,'uint16');
u_vec=zeros(3,3,max_iter,'uint16');
v_vec=zeros(3,3,max_iter,'uint16');

for agent_type_ind = 1:3
    agent_type = agent_types(agent_type_ind);
    parfor iter=1:max_iter
        % starting bids vector
        bid=zeros(1,3);
        
        % starting minimum prices (asks)
        min_p=ones(1,3)*800;

        % starting without matched/sold goods
        current_match = [0,0,0];

        % iterate over periods
        for i=1:T
            % draw the deviating player
            deviator_is_buyer=randi(2)-1; % =1 buyer, =0 seller
            ind=randi(3); % player index 1,2,or 3

            if deviator_is_buyer
                % first check if the player has a deviation. Two
                % conditions for buyers:
                % 1. Matched players cannot deviate
                if current_match(ind)==0
                    % all goods that give positive payoff to the deviator
                    if strcmp(agent_type,"uniform")
                        % for uniform agents any deviation to a
                        % positive payoff is possible
                        better_goods=find(h(:,ind)-max(min_p(:),bid(:))-1>=0);
                    else%
                        % in both other cases only the best goods (highest margin)
                        % are considered as possible deviations
                        max_val=max(h(:,ind)-max(min_p(:),bid(:))-1);
                        better_goods=find(h(:,ind)-max(min_p(:),bid(:))-1==max_val);
                    end
                    if ~isempty(better_goods) % 2. There is a better bid
                        if max(h(:,ind)-max(min_p(:),bid(:))-1)>=0
                            % pick a random good to bid for out of
                            % better_goods
                            g_ind=randi(length(better_goods));

                            % unmatch however is currently the highest
                            % bidder
                            current_match(current_match==better_goods(g_ind))=0;

                            if strcmp(agent_type,"best reply")
                                % increment the bid
                                bid(better_goods(g_ind))=max(min_p(better_goods(g_ind)),bid(better_goods(g_ind)))+1;
                            else
                                % in both other cases
                                % randomly draw the bid
                                if h(better_goods(g_ind),ind) - max(min_p(better_goods(g_ind)),bid(better_goods(g_ind)))>0
                                    bid(better_goods(g_ind))=randi([max(min_p(better_goods(g_ind)),bid(better_goods(g_ind)))+1, h(better_goods(g_ind),ind)]);
                                end
                            end
                            % record the buyer as matched to the good
                            current_match(ind)=better_goods(g_ind);
                        end
                    end
                end
            else
                % deviator is seller
                % check that the good is currently not sold
                if bid(ind)==0
                    % reduce the price without going below the cost
                    if min_p(ind)>c(ind)
                        min_p(ind)=randi([c(ind),min_p(ind)-1]);
                    end
                end
            end

            % check if the iteration has converged
            not_converged=false;

            % The iteration has not converged if:
            % 1. unmatched buyers have a deviation
            for b=1:3
                if current_match(b)==0
                    better_goods=find(h(:,b)-max(min_p(:),bid(:))-1>=0, 1);
                    if ~isempty(better_goods) % buyers have no deviation
                        not_converged=true;
                    end
                end
            end

            if ~all(min_p==c | bid>0) % 2. if sellers can reduce prices
                not_converged=true;
            end


            if ~not_converged
                % record final utilities
                u=zeros(3,1);
                v=zeros(3,1);
                for b=1:3
                    if current_match(b)>0
                        u(b)=h(current_match(b),b)-bid(current_match(b));
                        v(current_match(b))=bid(current_match(b))-c(current_match(b));
                    end
                end
                
                u_vec(agent_type_ind,:,iter) = u;
                v_vec(agent_type_ind,:,iter) = v;
                
                % match the third buyer and third seller for 0 payoff to
                % aggregate matchings like [210] and [213]
                if all(current_match == [1 2 0]) || all(current_match == [2 1 0])
                    current_match(3) = 3;
                end
                
                % record final match
                final_match(agent_type_ind,:,iter)=current_match;

                % record final prices
                final_prices(agent_type_ind,:,iter) = max(bid,min_p);
                
                break
            else
                if i == T % Maximum periods reached without convergence
                    disp('Error: simulation did not converge - consider increasing the number of rounds T');
                end
            end
        end
    end
end