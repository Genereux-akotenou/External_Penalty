// ****************************************************
// FUNCTION
// ****************************************************

function [val, grad] = f(x1, x2)
    val  = 2*(x1^2) + 3*(x1*x2) + 2*(x2^2)
    grad = [4*x1 + 3*x2; 4*x2 + 3*x1]
endfunction

function [penalty]=P(x)
    x1 = x(1)
    x2 = x(2)
    penalty = (max(0, x1 + 0.5))^2 + (max(0, x2 + 0.5))^2;
endfunction

function [val, grad] = phi(x, p)
    x1 = x(1); 
    x2 = x(2);
    [fx, grad_f] =  f(x1, x2);
    val  = fx + 0.5*p*P(x);
    grad = grad_f + [p*max(0, x1 + 0.5); p*max(0, x2 + 0.5)];
endfunction

function [Hv] = Hv(x, v, p)
    x1 = x(1); 
    x2 = x(2);
    H = zeros(2, 2)
    H(1, 1) = 4 + p*(x1 > -0.5)
    H(1, 2) = 3
    H(2, 1) = 3 
    H(2, 2) = 4 + p*(x2 > -0.5)
    Hv = v'*H;
endfunction

// ****************************************************
// Algorithme de la penalite exterieure
// ****************************************************

function [x_k, X_trace, P_trace, K_trace, J_trace] = ExternalPenalty(x0, p0, c, M, eps1, eps2, direction_type, MaxIter, verbose)
    // Init
    J_trace = [];
    K_trace = [];
    X_trace = [];
    P_trace = [];
    k   = 0
    x_k = x0
    p_k = p0
    
    // Armijo {0<τ0<0.5}
    tau = 0.25
    
    // For Newton_Direction_Only
    Delta=1.20
    
    // Let's display infos
    pad = '';
    if verbose>0,
        s='   ';
        for i=1:verbose,
          pad = pad + s
        end
        mprintf("\n%s DIRECTION_TYPE=[%s], c=%d, p0=%d\n%s --------------\n", pad, direction_type, c, p0, pad)
        mprintf("%s iter      j       p_k           x_k                        ||nabla φ(x_k, p_k)||\n",pad)
        [_, nabla_phi] = phi(x_k, p_k)
        mprintf("%s %3d       %s   %10.1f      [%f, %f]        %10.7f\n", pad, k, "_", p_k, x_k(1), x_k(2), norm(nabla_phi))
    end
    
    // Main loops
    while (p_k < M) & (abs(P(x_k)) > eps1)  & (k < MaxIter)
        j = 0
        x_j = x_k
        [val_phi, nabla_phi] = phi(x_j, p_k)
        
        // Second loops
        while norm(nabla_phi) > eps2
            [val_phi, nabla_phi] = phi(x_j, p_k)
            
            // ... Descent direction ...
            if direction_type == "-nabla_phi"
                d_j = -nabla_phi
            else
                d_j = GC_TR(nabla_phi', Hv, x_k, p_k, eps1, Delta, MaxIter, 0)
            end
            // ... Feasible step     ...
            theta = 1
            while phi(x_j+theta*d_j, p_k) - val_phi > (theta*tau*nabla_phi*d_j')
                theta = theta / 2
            end
            // ... Update x_j and j  ...
            x_j = x_j + theta*d_j
            j = j + 1
        end
        
        // Update x_k, p_k, k(iter)
        x_k = x_j
        p_k = c * p_k
        k = k + 1
        
        // Trace
        if verbose>0,
            mprintf("%s %3d     %3d     %10.1f      [%f, %f]       %10.7f\n", pad, k, j, p_k, x_k(1), x_k(2), norm(nabla_phi))
        end
        
        // Record to make plot later
        J_trace = [J_trace, j]
        K_trace = [K_trace, k]
        X_trace = [X_trace, x_k]
        P_trace = [P_trace, p_k]
    end
endfunction
