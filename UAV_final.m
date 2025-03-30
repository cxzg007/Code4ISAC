close all
clear all

% 参数设置
L = 3; % UAV 的数量
K = 3;% GU 的数量
N = 1; % 目标的数量
H = 30; % UAV 高度
Pmax = 43; % 传输功率 (dBm)
rho_0 = -30; % 参考距离的信道功率增益 (dB)
%Gamma = 0;  SINR 阈值 (dB)
sigma_0_squared = -60; % 通信和感知噪声功率 (dBm)
d_min = 10; % 最小 UAV 间距离 (m)
u_t_n = [-10, 2]; % 目标的坐标
u_c_k = [-32,52;-92,-88;86,-47]; % GU的坐标

% 单位换算
Pmax = (10^(Pmax / 10)) * 10^(-3); % 转换为 W
sigma_0_squared = (10^(sigma_0_squared / 10)) * 10^(-3); % 转换为 W
Gamma = 1; % SINR 阈值 
rho_0 = 10^(rho_0 / 10); %参考距离的信道功率增益 

% 初始化 UAV 位置和功率
u_r_l = [-32,52;-92,-88;86,-47]; % UAV 初始位置
p_l_opt = ones(L, 1) * Pmax / 2; % UAV 初始功率（设为一半的最大功率）
% delta_l_k = zeros(L, K); % 假设初始 delta_l_k 值
tolerance = 1e-7; % 设定收敛阈值
max_iter = 10; % 最大迭代次数
iter = 0; % 迭代计数
zeta_tmp = 0;

% 交替优化过程
while true
    iter = iter + 1;
    fprintf('Iteration: %d\n', iter);
    
    % 步骤 1：固定 UAV 位置，优化传输功率
    cvx_begin 
        variable p_l(L) nonnegative; 
        variable zeta;
        expression interference_sum;
        % 目标函数
        maximize(zeta)
       
        subject to
        % 约束条件 (6b)
            0 <= p_l <= Pmax;
            
            % 约束条件 (6c)
            for k = 1:K
                interference_sum = 0;
                for l = 1:L
                    if l ~= k
                        dist_l = norm(u_r_l(l,:) - u_c_k(k,:))^2 + H^2;
                        interference_sum = interference_sum + (p_l(l) * rho_0) / dist_l;
                    end
                end
                dist_s = norm(u_r_l(k,:) - u_c_k(k,:))^2 + H^2;
                p_l(k) * rho_0 / dist_s >= Gamma * (interference_sum + sigma_0_squared);
            end
            
            % 约束条件 (7a)
            for n = 1:N
                sensing_sum = 0;
                for l = 1:L
                    dist_sensing = sum_square(u_r_l(l,:) - u_t_n(n,:)) + H^2;
                    d_4 = pow_p(dist_sensing, 2);
                    sensing_sum = sensing_sum + p_l(l) / d_4;
                end
                sensing_sum >= zeta;
            end
    cvx_end

    % 保存功率优化结果
    p_l_opt = p_l;
    fprintf('p_l_opt:\n');
    disp(p_l_opt);

    % 步骤 2：固定传输功率，优化 UAV 位置
    cvx_begin 
        variables u_l(L, 2) delta_l_k(L,K) zeta;
        expression dist_sq
        expression f_ul
        maximize(zeta)

        subject to
            % 约束条件 (20a)
            for n = 1:N
                sum_terms = 0;
                for l = 1:L
                    dist_squared = sum_square(u_r_l(l,:) - u_t_n(n,:)) + H^2;
                    d_tilde_4 = pow_p(dist_squared, 2);
                    
                    
                    dist_sq = norm([u_l(l,:) - u_t_n(n,:),H]);
                    d_t_4 = pow_pos(dist_sq, 4);
                    
                    d_tilde_8 = pow_p(dist_squared, 4);
                    sum_terms = sum_terms + (2 * p_l_opt(l) / d_tilde_4 - (p_l_opt(l) * d_t_4) / d_tilde_8);
                end
                sum_terms >= zeta;
            end
            
            % 约束条件 (20b)
            for l = 1:L
                for i = 1:L
                    if l ~= i
                        norm(u_r_l(l, :) - u_r_l(i, :))^2 + 2 * ((u_r_l(l, :) - u_r_l(i, :))' * (u_l(l, :) - u_l(i, :))) >= d_min^2;
                    end
                end
            end
            
            % 约束条件 (20c)
            for l = 1:L
                sum_term = 0;
                for i = 1:L
                    if l ~= i
                        sum_term = sum_term + (p_l_opt(i) * rho_0) * inv_pos(H^2 + delta_l_k(i,k));
                    end
                end
                sum_term = sum_term + sigma_0_squared;
            end
           for k = 1:K
                dist = sum_square(u_r_l(k, :) - u_c_k(k, :)) + H^2;
                f_ul = p_l_opt(k) * rho_0 / dist - (p_l_opt(k) * rho_0) / (dist^2) * (sum_square(u_l(k,:) - u_c_k(k,:)) - sum_square(u_r_l(k,:) - u_c_k(k,:)));
                f_ul - Gamma * sum_term >= 0;
            end
            
            % 约束条件 (20d)
            for l = 1:L
                for k = 1:K
                    if l ~= k
                        norm(u_r_l(l, :) - u_c_k(k, :))^2 + 2 * ((u_r_l(l, :) - u_c_k(k, :))' * (u_l(l, :) - u_r_l(k, :))) >= delta_l_k(l, k);
                    end
                end
            end
    cvx_end
    
    % 更新 UAV 位置
    u_r_l = u_l;
    fprintf('u_r_l\n');
    disp(u_r_l);

    % 检查收敛性
    if (zeta - zeta_tmp)< tolerance || iter >= max_iter
        break;
        
    zeta_tmp=zeta;
    
    end
end

% 输出结果
fprintf('Optimized UAV positions:\n');
disp(u_r_l);
fprintf('Optimized UAV transmit powers:\n');
disp(p_l_opt);




