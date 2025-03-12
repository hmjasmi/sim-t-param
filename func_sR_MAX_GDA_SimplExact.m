function [sum_rate, opt] = func_sR_MAX_GDA_SimplExact(H_RI, H_IT, Tl, N_0, properties)
% This function maximizes the sum rate of the system using GDA
% The channel model for optimization is Simplified
% The sum_rate returned for the optimal solution considers Simplified channel for computation
% The opt.sum_rate_simpleExact returned for the optimal solution considers Exact channel for computation


% The channel model for returned optimal channel is Exact 
    maxIter = properties.GDA.maxIter;
    tol     = properties.GDA.tol;
    alpha   = properties.GDA.alpha;  % initial step size
    delta   = properties.GDA.delta;
    tau     = properties.GDA.tau;
    K       = properties.No_Users;
    N       = properties.No_SIM_RE_z * properties.No_SIM_RE_y;
    L       = properties.No_SIM_layers;
    Sl_21   = properties.Sl_blocks.S21;
    % Initialization
    [Theta, G_matrices, debug] = func_MRT_init(H_RI, H_IT, properties);
    theta = [];
    for l = 1:L
        theta = [theta; angle(diag(Theta{l}).')];
    end

    iter = 1;
    while true
        % Compute current sum rate based on theta (using the helper function)
        Tl_simpl = zeros(2*N,2*N);
        Tl_simpl(N+1:end,N+1:end) = inv(Sl_21);
        f_old = compute_sum_rate(theta, H_RI, H_IT, Tl_simpl, N, L, N_0);
        sum_rate(iter) = f_old;
        sum_rate_simpleExact(iter) = compute_sum_rate(theta, H_RI, H_IT, Tl, N, L, N_0);
        % Compute numerical gradient
        grad = zeros(L, N);
        for l = 1:L
            for n = 1:N
                theta_plus = theta;
                theta_plus(l, n) = theta_plus(l, n) + delta;
                f_plus = compute_sum_rate(theta_plus, H_RI, H_IT, Tl_simpl, N, L, N_0);
                grad(l, n) = -(f_plus - f_old) / delta;
            end
        end
        
        
        % --- Armijo Backtracking Line Search ---
        c = 1e-4;  % sufficient increase constant
        theta_old = theta;  % store current theta
        
        % Compute candidate update with current alpha (do not update theta yet)
        theta_cand = theta_old - alpha * grad;
        f_cand = compute_sum_rate(theta_cand, H_RI, H_IT, Tl_simpl, N, L, N_0);
        
        % Backtracking: reduce alpha until the Armijo condition is met
        while f_cand < f_old + c * alpha * norm(grad, 'fro')^2
            alpha = tau * alpha;  % tau < 1 reduces alpha
            theta_cand = theta_old - alpha * grad;
            f_cand = compute_sum_rate(theta_cand, H_RI, H_IT, Tl_simpl, N, L, N_0);
            if alpha < 1e-8
                fprintf('Alpha became very small during backtracking.\n');
                break;
            end
        end
        
        % Update theta once a satisfactory step is found
        theta = theta_cand;
        fprintf('Iteration %d, sR: %e, Alpha: %e\n', iter, f_old, alpha);
        
        % Optionally, reset alpha if desired (or let it carry over)
        alpha = properties.GDA.alpha;  % reset to initial value
        
        % Check convergence: relative improvement or max iterations
        if iter > 1 && abs(f_cand - f_old)/abs(f_old) < tol || (iter > maxIter)
            break;
        end
        
        iter = iter + 1;
    end

    % Store optimal values
    theta_opt = theta;
    for l = 1:L
        G_matrices_opt{l} = blkdiag(diag(exp(1j * theta_opt(l,:))), diag(exp(-1j * theta_opt(l,:))));
        Theta_opt{l} = diag(exp(1j * theta_opt(l,:)));
    end

    T_I_opt = G_matrices_opt{1} * Tl;
    for l = 2:L-1
        T_I_opt = T_I_opt * G_matrices_opt{l} * Tl;
    end
    T_I_opt = T_I_opt * G_matrices_opt{L};
    T_I_22_opt = T_I_opt(N+1:end, N+1:end);

    T_I_opt_ss = G_matrices_opt{1} * Tl_simpl;
    for l = 2:L-1
        T_I_opt_ss = T_I_opt_ss * G_matrices_opt{l} * Tl_simpl;
    end
    T_I_opt_ss = T_I_opt_ss * G_matrices_opt{L};
    T_I_22_opt_ss = T_I_opt_ss(N+1:end, N+1:end);
    H_opt_ss = H_RI * inv(T_I_22_opt_ss) * H_IT;
    H_opt = H_RI * inv(T_I_22_opt) * H_IT;

    opt.Theta_opt = Theta_opt;
    opt.G_matrices_opt = G_matrices_opt;
    opt.H_opt = H_opt;
    opt.H_opt_ss = H_opt_ss;
    opt.T_I_opt = T_I_opt;
    opt.T_I_opt_ss = T_I_opt_ss;
    opt.T_I_22_opt = T_I_22_opt;
    opt.T_I_22_opt_ss = T_I_22_opt_ss;
    opt.sum_rate_simpleExact = sum_rate_simpleExact;
end

function f = compute_sum_rate(theta, H_RI, H_IT, Tl, N, L, N_0)
    % Build the G_matrices based on theta
    for ll = 1:L
        G_matrices{ll} = blkdiag(diag(exp(1j * theta(ll,:))), diag(exp(-1j * theta(ll,:))));
    end
    T_I = G_matrices{1} * Tl;
    for ll = 2:L-1
        T_I = T_I * G_matrices{ll} * Tl;
    end
    T_I = T_I * G_matrices{L};
    T_I_22 = T_I(N+1:end, N+1:end);
    H = H_RI * inv(T_I_22) * H_IT;
    [P_R, P_I] = func_signalPower(abs(H));
    SINR = P_R ./ (sum(P_I,2) + N_0);
    f = sum(log2(1 + SINR), 1);
end
