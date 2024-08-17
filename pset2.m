% Define the quantized state and control levels
states = [0.0, 0.5, 1.0];
controls = [-0.4, -0.2, 0.0, 0.2, 0.4];

% Initialize the cost-to-go function J2
J2 = 4 * abs(states);

% Initialize the cost-to-go function J1 and optimal control for k=1
J1 = zeros(1, length(states));
u1_opt = zeros(1, length(states));

% Compute J1 and u1_opt
for i = 1:length(states)
    x1 = states(i);
    min_cost = inf;
    
    for u1 = controls
        x2 = x1 - 0.4 * x1^2 + u1;
        
        % Interpolate J2 for non-quantized x2 values
        if x2 < 0
            continue; % Skip out-of-bounds x2 values
        elseif x2 > 1.0
            continue; % Skip out-of-bounds x2 values
        else
            J2_interp = interp1(states, J2, x2, 'linear');
            % disp("X2 : ");
            % disp(x2);
            % disp("J2 : ");
            % disp(J2_interp);
        end
        
        cost = abs(u1) + J2_interp;
        
        %disp("Cost : ");
        % disp(cost);
        
        if cost < min_cost
            min_cost = cost;
            u1_opt(i) = u1;
        end
    end
    
    J1(i) = min_cost;
end

% Display the results for k=1
fprintf('Optimal control u1 and cost J1 for each state x1:\n');
for i = 1:length(states)
    fprintf('x1 = %.1f, u1 = %.1f, J1 = %.2f\n', states(i), u1_opt(i), J1(i));
end

% Initialize the cost-to-go function J0 and optimal control for k=0
J0 = zeros(1, length(states));
u0_opt = zeros(1, length(states));

% Compute J0 and u0_opt
for i = 1:length(states)
    x0 = states(i);
    min_cost = inf;
    
    for u0 = controls
        x1 = x0 - 0.4 * x0^2 + u0;
        
        % Interpolate J1 for non-quantized x1 values
        if x1 < 0
            continue; % Skip out-of-bounds x1 values
        elseif x1 > 1.0
            continue; % Skip out-of-bounds x1 values
        else
            J1_interp = interp1(states, J1, x1, 'linear');
            disp("X1 : ");
            disp(x1);
            disp("J1 : ");
            disp(J1_interp);        
        end
        cost = abs(u0) + J1_interp;
        
        if cost < min_cost
            min_cost = cost;
            u0_opt(i) = u0;
        end
    end
    
    J0(i) = min_cost;
end

% Display the results for k=0
fprintf('Optimal control u0 and cost J0 for each state x0:\n');
for i = 1:length(states)
    fprintf('x0 = %.1f, u0 = %.1f, J0 = %.2f\n', states(i), u0_opt(i), J0(i));
end


% Parameters
A = 1.05;    % System matrix
B = 1;       % Input matrix
q = 1;       % Weight on state (adjust as needed)
r = 1;       % Weight on control (adjust as needed)

% Calculate the steady-state solution for P using the DARE
P = dare(A, B, q, r);

% Calculate the optimal gain F
F = (A * P) / (r + B * P);

% Display the results
fprintf('Steady-State Solution P: %.4f\n', P);
fprintf('Optimal Gain F: %.4f\n', F);

A = [1 1; 1 0];
B = [1; 0];
Q = [1 0; 0 1];
R = 1;

% Solve the Discrete Algebraic Riccati Equation
P_ss = dare(A, B, Q, R);

disp('Steady-State Solution P_ss:');
disp(P_ss);
