function dist=find_distance_to_lcore(x_point,solv,dist_mode)
% use the optimizer with the provided point
[res, prob]=solv(x_point);

if prob~=0 % Problem in the solver
    disp('Error: could not solve the optimization problem');
end

% return the distance between x_point and the closest point in the L-core
dist=find_distance_to_point(double(x_point),res{3},dist_mode);