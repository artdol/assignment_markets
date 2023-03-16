function dist=find_distance_to_point(x_point,p,dist_mode)
% finds distance between x_point and p using one of the three distance
% regimes: "euc" (euclidean), "chess" or "city" (Manhattan or cityblock)

if strcmp(dist_mode,"euc")
% Euclidian
dist=norm(x_point-p);

elseif strcmp(dist_mode,"chess")
% Chessboard chebyshev
dist=max(abs(x_point-p));

elseif strcmp(dist_mode,"city")
% Manhattan cityblock
dist=sum(abs(x_point-p));
end
end