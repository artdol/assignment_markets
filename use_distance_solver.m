function dist = use_distance_solver(solv,x_point,dist_mode)
% solv = optimizer object
% x_point = the point to which the distance is calculated
% dist_mode = distance regime

% solve the optimization program
[sol, res]=solv(x_point);

if res~=0
    disp('Error: could not solve the optimization problem');
end

% return the distance between x_point and the closest point in the core
dist=find_distance_to_point(double(x_point),sol{3},dist_mode);
