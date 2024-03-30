% Simulate the Pairwise Nash core for random weights
% This produces the table of frequencies of matchings (Table 5)

% number of simulated data points
max_iter = 50000;

% Small positive epsilon for strict inequalities
eps=0.0001;

% generate all full matchings
all_matchings=perms(1:3);
num_matchings=length(all_matchings);

% generate data points and check, which matchings are stable. If there are
% multiple matchings that are solutions for the same weights then count each as an individual point.
is_stable=zeros(max_iter,num_matchings);
for iter=1:max_iter
    % randomly draw weights
    alphas=rand(1,5);
    alphas=diff(sort([0 alphas 1])); % this ensures that the draw is uniform
    
    % iterate over full matchings
    for i=1:num_matchings
        m=all_matchings(i,:);
        
        % calculate the utilities given the matching and weights
        buyer_u=zeros(1,3);
        seller_u=zeros(1,3);
        for b=1:3
            buyer_u(b)=A(m(b),b)*alphas(3+b)/(alphas(3+b)+alphas(m(b)));
            seller_u(m(b))=A(m(b),b)*alphas(m(b))/(alphas(3+b)+alphas(m(b)));
        end
        
        % check if players can deviate
        is_dev=false;
        for b=1:3
            for s=1:3
                if s~=m(b)
                    buyer_p=A(s,b)*alphas(3+b)/(alphas(3+b)+alphas(s));
                    seller_p=A(s,b)*alphas(s)/(alphas(3+b)+alphas(s));

                    if buyer_p>=buyer_u(b) && seller_p>=seller_u(s)
                        if buyer_p>=buyer_u(b) +eps || seller_p>=seller_u(s) +eps
                            is_dev=true;
                        end
                    end
                end
            end
        end
        if ~is_dev
            % if not, record matching as stable, i.e. Pairwise Nash core
            % solution
            is_stable(iter,i)=1;
        end
    end
end

% relative frequency of stable matchings
res=sum(is_stable)/max_iter;

% sort by frequency
[sorted_pairings_cnt, pairings_order] = sortrows(res','descend');
sorted_pairings=all_matchings(pairings_order,:);

% Print the table of matchings and frequencies
% Heading
fprintf('%s & %s \\\\ \n',"Matching", "Frequency")
for i=1:num_matchings
    fprintf('$\\mathit{[%d %d %d]}$ & %0.2f\\%% \\\\ \n',sorted_pairings(i,:),sorted_pairings_cnt(i)*100)
end
