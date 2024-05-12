// ***************************************************************************************
// UTILS FOR DESCENT DIRECTION
// ***************************************************************************************

function [disc, theta_max] = max_deplacement(d, delt, d_confiance)
    // -----------------------------------------------------------------------------------
    // JOB: calcul le pas de d́eplacement maximum pour demeurer dans la ŕegion de confiance
    // -----------------------------------------------------------------------------------
    A = d'*d
    B = 2*delt'*d
    C = delt'*delt-(d_confiance^2)
    disc = B^2 - 4*A*C
    theta_max = -1
    if disc >= 0 then
        theta_max = (-B + sqrt(disc)) / (2*A)
        if theta_max < 0 then
            theta_max = 1
        end
    end
endfunction

function [dN, dNQ, niter, L_f] = GC_TR(c, Hv, x, rho, eps, Delta, MaxIter, verbose)
    // -----------------------------------------------------------------------------------
    // JOB: Gradient conjugué_Région de confiance
    // -----------------------------------------------------------------------------------
    
    // init
    niter = 0; 
    L_f = [];
    [bidon,dim] = size(c)
    n = length(x)
    delt = zeros(n,1);
    deltQ = Hv(x, delt, rho);
    nablaq = c + deltQ;
  
    // Let's display infos
    pad = '';
    if verbose >0,
        s='   ';
        for i=1:verbose,
            pad = pad + s
        end
        mprintf("%sDEBUG=[GC_TR]", pad)
        mprintf("%s iter    q(delt)  ||nabla q(delt)||  theta\n",pad)
        mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
            norm(nablaq), 0)
    end
  
    // Let's compute ||∇q(x)||
    norm2nablaq = nablaq*nablaq';
    pasprecis = norm2nablaq>eps^2;
    sortie = %F
  
    // Let's compute d_0
    d = -nablaq'
    b = 0
    while pasprecis & (niter < MaxIter) & ~sortie
        // d_k
        niter = niter + 1;
        d = -nablaq' + (b*d);
        dQ = Hv(x, d, rho);
        dQd = dQ*d;
        
        // theta_max & courbure_negative
        [disc, theta_max] = max_deplacement(d, delt, Delta)
        sortie = (dQd <= 0) || (disc<0)
    
        if ~sortie then
            // theta_k
            theta=(-nablaq*d)/dQd;
            if theta<0 | theta>theta_max then
                theta = theta_max
            end
            
            // delt_k+1
            delt = delt + theta*d;
            
            // update nabla_q, deltQ, beta
            nablaq = nablaq + theta*dQ;
            deltQ = deltQ + theta*dQ;
            
            nablaq_T = nablaq'
            nablaq_Q = Hv(x, nablaq_T, rho)
            b = nablaq_Q*d / dQd
           // b = dQ*nablaq'/ dQd
            
            // update stop conditions    
            norm2nablaq = nablaq*nablaq';
            pasprecis = norm2nablaq>eps^2;
            
            if verbose>0,
            mprintf("%s %3d  %10.7f    %10.7e   %10.7f\n", pad, niter, (c*delt + 0.5*deltQ*delt), ...
                norm(nablaq), theta)
            end
            L_f(niter) = c*delt + 0.5*deltQ*delt;
        end
    end
    dN = delt;
    dNQ = deltQ;
endfunction
