// ****************************************************
// WE RUN PENALTY ALGORITHM HERE
// ****************************************************

// Problem data
x0 = [-0.3; 0.5]
rho = 1
c  = 10
M  = 10^9
eps1 = 10^-6
eps2 = 10^-6
verbose = 1
MaxIter = 20

//Let's import our built in function
exec('PNE_search.sce', -1);
exec('PNE_Newton_direction.sce', -1);

// ****************************************************
// USE DIRECTION AS -nabla_phi FOR c=10 and 
// ****************************************************
direction_type = "-nabla_phi"
[x_opt, X_trace, P_trace, K_trace, J_trace] = ExternalPenalty(x0, rho, 10, M, eps1, eps2, direction_type, MaxIter, verbose)
[x_opt, X_trace, P_trace, K_trace, J_trace] = ExternalPenalty(x0, rho, 2, M, eps1, eps2, direction_type, MaxIter, verbose)

// ****************************************************
// LET USE DIRECTION_NEWTON_MODIFIÉE FOR c=10 and 
// ****************************************************
direction_type = "Direction_Newton_Modifiée"
[x_opt, X_trace, P_trace, K_trace, J_trace] = ExternalPenalty(x0, rho, 10, M, eps1, eps2, direction_type, MaxIter, verbose)
[x_opt, X_trace, P_trace, K_trace, J_trace] = ExternalPenalty(x0, rho, 2, M, eps1, eps2, direction_type, MaxIter, verbose)

// Plot
