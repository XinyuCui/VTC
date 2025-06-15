clc;clear all;clf;
close all;

%由于各时隙的卫星拓扑是固定的，如果每跑一次程序都要重新生成一遍的话，就会浪费时间
%所以该程序直接将各时隙的adjacency_matrix生成好，供后续程序直接调用

% 步骤 1: 参数设置
T = 66; % 卫星总数
P = 6; % 轨道平面数
F = 2; % 相位因子
N = T / P; % 每个轨道平面的卫星数

% 轨道参数
radius = 7000; % 轨道半径 (km)
inclination = 86.4; % 轨道倾角 (度)
inclination = deg2rad(inclination); % 转换为弧度

mu = 3.986004418e14; % 地球引力常数 (m^3/s^2)
n = sqrt(mu / (radius*1000)^3); % 平均运动角速度 (rad/s)

% 时间参数
total_time = 240; % 总模拟时间 (s)
time_step = 2; % 时间步长 (s)
num_steps = total_time / time_step;

% 用户坐标
user_coords = [6370, 0, 0; -6370, 0, 0];
num_users = size(user_coords, 1);

adjacency_matrix_of_each_time_step = false(num_steps, T + num_users, T + num_users);
sat_positions_of_each_time_step=zeros(num_steps, T + num_users, 3);

% 循环计算每个时间步的拓扑
for step = 1:num_steps
    t = (step - 1) * time_step;

    % 步骤 2: 卫星位置计算
    sat_positions = zeros(T, 3);
    for j = 1:P
        RAAN = 2 * pi * (j - 1) / P; % 升交点赤经
        for k = 1:N
            % 计算平近点角，考虑时间变化
            M = 2 * pi * (k - 1) / N + 2 * pi * F * (j - 1) / T + t * n;
            % 计算卫星在轨道平面内的位置
            x = radius * cos(M);
            y = radius * sin(M);
            z = 0;
            % 旋转到轨道平面
            R_RAAN = [cos(RAAN) -sin(RAAN) 0;
                      sin(RAAN) cos(RAAN) 0;
                      0 0 1];
            R_inclination = [1 0 0;
                             0 cos(inclination) -sin(inclination);
                             0 sin(inclination) cos(inclination)];
            pos = [x; y; z];
            pos = R_RAAN * R_inclination * pos;
            sat_index = (j - 1) * N + k;
            sat_positions(sat_index, :) = pos';
        end
    end

    % 步骤 3: 连接判断
    % 扩展邻接矩阵以包含用户
    adjacency_matrix = false(T + num_users, T + num_users);
    for j = 1:P % 轨道平面 index
        for k = 1:N % 每个轨道平面的卫星数 index
            current_sat = (j - 1) * N + k;

            % 前一个卫星
            prev_sat_in_orb = (j - 1) * N + mod(k - 2, N) + 1;
            adjacency_matrix(current_sat, prev_sat_in_orb) = true;
            adjacency_matrix(prev_sat_in_orb, current_sat) = true;

            % 后一个卫星
            next_sat_in_orb = (j - 1) * N + mod(k, N) + 1;
            adjacency_matrix(current_sat, next_sat_in_orb) = true;
            adjacency_matrix(next_sat_in_orb, current_sat) = true;

            % 左轨道卫星
            left_orb = mod(j - 2, P) + 1;
            left_sat_start = (left_orb - 1) * N + 1;
            left_sat_end = left_orb * N;
            left_sat_distances = zeros(N, 1);
            for left_k = 1:N
                left_sat = (left_orb - 1) * N + left_k;
                left_sat_distances(left_k) = norm(sat_positions(current_sat, :) - sat_positions(left_sat, :));
            end
            [~, nearest_left_index] = min(left_sat_distances);
            nearest_left_sat = (left_orb - 1) * N + nearest_left_index;
            adjacency_matrix(current_sat, nearest_left_sat) = true;
            adjacency_matrix(nearest_left_sat, current_sat) = true;

            % 右轨道卫星
            right_orb = mod(j, P) + 1;
            right_sat_start = (right_orb - 1) * N + 1;
            right_sat_end = right_orb * N;
            right_sat_distances = zeros(N, 1);
            for right_k = 1:N
                right_sat = (right_orb - 1) * N + right_k;
                right_sat_distances(right_k) = norm(sat_positions(current_sat, :) - sat_positions(right_sat, :));
            end
            [~, nearest_right_index] = min(right_sat_distances);
            nearest_right_sat = (right_orb - 1) * N + nearest_right_index;
            adjacency_matrix(current_sat, nearest_right_sat) = true;
            adjacency_matrix(nearest_right_sat, current_sat) = true;
        end
    end

    % 检查用户与卫星的连接
    for u = 1:num_users
        user_pos = user_coords(u, :);
        distances = vecnorm(sat_positions - repmat(user_pos, T, 1), 2, 2);
        close_sats = find(distances <= 2000); %用户可与相距2000km的卫星保持连接
        for sat = close_sats
            user_index = T + u;
            adjacency_matrix(user_index, sat) = true;
            adjacency_matrix(sat, user_index) = true;
        end
    end

    % 步骤 4: 图的构建与更新
    G = graph(adjacency_matrix);

%     % 步骤 5: 3D 可视化
%     figure;
%     % 绘制卫星节点
%     scatter3(sat_positions(:,1), sat_positions(:,2), sat_positions(:,3), 'b', 'filled');
%     hold on;
%     % 绘制连接边
%     [s, t] = find(adjacency_matrix(1:T, 1:T));
%     for i = 1:length(s)
%         start_point = sat_positions(s(i), :);
%         end_point = sat_positions(t(i), :);
%         line([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'Color', 'r');
%     end
%     % 绘制用户节点
%     scatter3(user_coords(:,1), user_coords(:,2), user_coords(:,3), 'g', 'filled');
%     % 绘制用户与卫星的连接
%     for u = 1:num_users
%         user_index = T + u;
%         connected_sats = find(adjacency_matrix(user_index, 1:T));
%         for sat = connected_sats
%             user_pos = user_coords(u, :);
%             sat_pos = sat_positions(sat, :);
%             line([user_pos(1), sat_pos(1)], [user_pos(2), sat_pos(2)], [user_pos(3), sat_pos(3)], 'Color', 'g');
%         end
%     end
%     hold off;
%     title(sprintf('Walker Delta 星座拓扑结构 ( t = %d s 至 t = %d s )', (step-1)*time_step, step*time_step));
%     xlabel('X (km)');
%     ylabel('Y (km)');
%     zlabel('Z (km)');
%     grid on;

    % 将 user_coords 的第一行加入到 sat_positions 的倒数第二行，将 user_coords 的最后一行加入到 sat_positions 的最后一行
    % 这样sat_positions的倒数第二行就是远端视频库的位置，最后一行就是用户的位置了
    first_user_row = user_coords(1, :);
    last_user_row = user_coords(end, :);
    sat_positions = [sat_positions; first_user_row; last_user_row];
        
    %步骤6：将各时隙的adjacency_matrix存储到adjacency_matrix_of_each_time_step中
    %将各时隙的sat_positions存储到sat_positions_of_each_time_step中
    adjacency_matrix_of_each_time_step(step,:,:)=adjacency_matrix;
    sat_positions_of_each_time_step(step,:,:)=sat_positions;
end

%adjacency_matrix_of_each_time_step(1,:,:)是一个（1，68，68）维的数组
%在用的时候注意用squeeze函数对维度压缩一下，压缩成（68,68）维的，squeeze函数的作用就是去除数组中维度为1的部分
%同理sat_positions也是
save('adjacency_matrix_of_each_time_step.mat', 'adjacency_matrix_of_each_time_step');
save('sat_positions_of_each_time_step.mat', 'sat_positions_of_each_time_step');