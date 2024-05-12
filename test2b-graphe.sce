// ****************************************************
// FUNCTION
// ****************************************************

// Return f(x) and ∇f(x)  -> [scalar, vect(2)]
function [val, grad] = f(x1, x2)
    val  = 2*(x1^2) + 3*(x1*x2) + 2*(x2^2)
    grad = [2*x1 + 3*x2; 2*x2 + 3*x1]
endfunction

// Return φ(x1, x2, ρ) -> [scalar, vect(2)]
function [val, grad] = phi_value(x1, x2, p)
    [fx, grad_f] =  f(x1, x2);
    val  = fx + 0.5*p*(max(0, x1 + 0.5))^2 + 0.5*p*(max(x2 + 0.5))^2;
    grad = grad_f + [p*max(0, x1 + 0.5); p*max(0, x2 + 0.5)];
endfunction

// Return surface φ(., ρ) -> matrix(x1.length, x2.length)
function [val] = phi_surface(x1_values, x2_values, p)
    n1 = length(x1_values)
    n2 = length(x2_values)
    
    val  = zeros(n1, n2)
    for i = 1:n1
        x1 = x1_values(i)
        for j = 1:n2
            x2 = x2_values(j)
            [fx, grad_f] =  f(x1, x2)
            val(i, j)  = fx + 0.5*p*(x1 + 0.5)^2 + 0.5*p*(x2 + 0.5)^2
        end
    end
endfunction

// Return ∇φ(p) -> vector[x1, x2]
function [grad] = nabla_phi(p)
    grad = zeros(2)
    x2 = (p^2/2 + (5*p)/2) / (9 - (2 + p)^2) 
    x1 = (-3*x2 - p/2) / (2 + p)
    grad(1) = x1
    grad(2) = x2
endfunction

// ****************************************************
// PLOT
// ****************************************************

// Define range for x1 and x2
x1_range = linspace(-1.5, 1.5, 10);
x2_range = linspace(-1.5, 1.5, 10);

// Define values of rho
rho_values = [5, 10, 100, 1000];

// Gradient points for each value of rho
gradient_points = [];
for rho = rho_values
    grad = nabla_phi(rho);
    disp(grad)
    gradient_points = [gradient_points, grad];
end

// Plot contour of phi & gradient points
//clf
//scf(0)
figure('position', [100, 100, 1000, 800]);
for i = 1:length(rho_values)
    rho = rho_values(i)
    surface = phi_surface(x1_range, x2_range, rho)
    contour(x1_range, x2_range, surface, 7, 'label', string(rho))
    scatter(gradient_points(1,:), gradient_points(2,:), '+');
end
xlabel('x1');
ylabel('x2');
title('Contours de la fonction de coût phi(., rho) pour différentes valeurs de rho');
