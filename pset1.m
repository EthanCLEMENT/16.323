fun = @(x) 2.*x(1).^4 + 3.*x(2).^3 + 6.*x(3).^2 -3.*x(1)*x(2) -6.*x(2)*x(3);
val=fminsearch(fun,[0 0 0]);
min=fun(val);

disp("F min : ");
disp(min);

syms x y z
f = 2*x^4 + 3*y^2 + 6*z^2 - 3*x*y - 6*y*z;
H = hessian(f,[x,y,z]);


disp("Hessian matrix : ");
disp(H);

Heig = [24*0.62^2 -3 0; -3 6 -6; 0 -6 12];

disp("Eigenvalues of H : ");
disp(eig(Heig));

syms x y z
f = 2*x^4 + 3*y^3 + 6*z^2 - 3*x*y - 6*y*z;

% Calculate the gradient symbolically
grad_f = gradient(f, [x, y, z]);

% Convert symbolic expressions to MATLAB functions
grad_f_fun = matlabFunction(grad_f, 'Vars', [x, y, z]);
f_fun = matlabFunction(f, 'Vars', [x, y, z]);

x0 = [1; 1; 1]; % Initial guess
tolerance = 1e-6; % Convergence tolerance
iter_max = 10000; % Maximum number of iterations
alpha_values = [0.1, 0.2, 0.5, 1]; % Different step sizes

% Loop over different alpha values
for alpha = alpha_values
    x = x0; % Reset initial value
    iter = 0; % Reset iteration counter
    grad = grad_f_fun(x(1), x(2), x(3)); % Initial gradient

    while (norm(grad, 2) >= tolerance)
        x_new = x - alpha * grad; % Update solution
        x = x_new; % Update old solution
        grad = grad_f_fun(x(1), x(2), x(3)); % Recompute the gradient
        iter = iter + 1;

        if iter > iter_max
            disp(['Maximum iterations reached for alpha = ', num2str(alpha)]);
            break;
        end
    end

    disp(['Alpha: ', num2str(alpha)]);
    disp(['Minimum found at: ', num2str(x')]);
    disp(['Function value at minimum: ', num2str(f_fun(x(1), x(2), x(3)))]);
    disp(['Number of iterations: ', num2str(iter)]);
    disp('-----------------------------');
end

syms x y lambda1 lambda2 lambda3

F = x^2 + y^2 - 6 * x * y - 4 * x - 5 * y;

f1 = -2 * x + y + 1;
f2 = x + y - 4;
f3 = x + 1;

L = F + lambda1 * f1 + lambda2 * f2 + lambda3 * f3;

dL_dx = diff(L, x);
dL_dy = diff(L, y);
dL_dlambda1 = diff(L, lambda1);
dL_dlambda2 = diff(L, lambda2);
dL_dlambda3 = diff(L, lambda3);

eqns = [dL_dx == 0, dL_dy == 0, f1 >= 0, f2 <= 0, f3 >= 0, ...
    lambda1 * f1 == 0, lambda2 * f2 == 0, lambda3 * f3 == 0];
sol = solve(eqns, [x, y, lambda1, lambda2, lambda3]);

x_sol = double(sol.x);
y_sol = double(sol.y);
lambda1_sol = double(sol.lambda1);
lambda2_sol = double(sol.lambda2);
lambda3_sol = double(sol.lambda3);

disp('Solution:')
disp("X : ");
disp(x_sol);
disp("Y : ");
disp(y_sol);
disp("Lambda1 : ");
disp(lambda1_sol);
disp("Lambda2 : ");
disp(lambda2_sol);
disp("Lambda3 : ");
disp(lambda3_sol);
