clear;
clc;

% Define the orignal and final rotation matrices of the endpoint
% (定义末端初始和最终旋转矩阵，即由旋转变换矩阵得到的姿态矩阵）
R0 = [0 1 0; 0 0 -1; -1 0 0]; % Orignal rotation matrix
Rt = [1 0 0; 0 0 -1; 0 1 0];  % Final rotation matrix

% Define the orignal and final position vectors of the endpoint
% (定义末端初始和最终位置向量）
P0 = [0.0 0.20 -0.4970]';
Pt = [0.0 0.20 -0.199]'; % [0,0,0,-pi/2,0,0,0]

% Convert rotation matrix to Quaternion (将旋转矩阵转换为四元数)
% Call Matlab Function 'rotm2quat'
q0 = rotm2quat(R0); % q0 and qt are 1*4 vectors
qt = rotm2quat(Rt);

% Number of interplations (插值个数）
n_insert = 100;                  
q_Quater = zeros(n_insert, 4); % Used to save interpolated Quaternions
R_Inter  = zeros(3, 3, n_insert);
R        = zeros(3, 3, n_insert+1);
axes     = zeros(n_insert+1, 3);
angs     = zeros(n_insert+1, 1);
q_Inter  = zeros(7,n_insert+1);
D_q      = zeros(7,n_insert+1);
D_ang    = zeros(n_insert+1, 3);
q        = [0,0,0,-pi/2,0,0,0]';

% Implement Slerp (Spherical Linear Interpolation) (实施Slerp插值)
for i = 1:n_insert
    s = i / (n_insert + 1);    % 插值n_insert个，则将均分为(n_insert+1)份
    theta = acos(dot(q0, qt)); % Angle between q0 and qt
    % Apply the Slerp Equation to achieve interpolation
    q_Quater(i, :) = (sin((1-s)*theta) / sin(theta)) * q0 + (sin(s*theta) / sin(theta)) * qt;
    
    % Convert Quaternion parameters to Rotation matrix (将插值的四元数转换为旋转变换矩阵)
    % Call Matlab Function 'quat2rotm'
    R_Inter(:, :, i) = quat2rotm(q_Quater(i, :));
end


% Calculate samll displacements and angles matrix

% Calculate samll displacements
D_x = (Pt(1) - P0(1))/(n_insert + 1);
D_y = (Pt(2) - P0(2))/(n_insert + 1);
D_z = (Pt(3) - P0(3))/(n_insert + 1);


for i = 1:n_insert+1

    % Calculate the axis-angle parameters of the intermediate transformation
    if i==1
        R(:, :, i) = R0'*R_Inter(:, :, i);                  % 相对于R_{0}坐标系，计算从R_{0}到第一个插值矩阵R_{1}的变换矩阵
    elseif i==n_insert+1
        R(:, :, i) = R_Inter(:, :, i-1)'*Rt;                % 相对于R_{n_insert}坐标系，计算从R_{n_insert}到目标旋转矩阵R_{t}的变换矩阵   
    else
        R(:, :, i) = R_Inter(:, :, i-1)'*R_Inter(:, :, i);  % 相对于R_{i-1}坐标系，计算从R_{i-1}到第一个插值矩阵R_{i}的变换矩阵
    end

     % Calculate the axis-angle parameters of the intermediate transformation
    % if i==1
    %     R(:, :, i) = R_Inter(:, :, i)*R0';                  % 相对于R_{0}坐标系，计算从R_{0}到第一个插值矩阵R_{1}的变换矩阵
    % elseif i==n_insert+1
    %     R(:, :, i) = Rt*R_Inter(:, :, i-1)';                % 相对于R_{n_insert}坐标系，计算从R_{n_insert}到目标旋转矩阵R_{t}的变换矩阵   
    % else
    %     R(:, :, i) = R_Inter(:, :, i)*R_Inter(:, :, i-1)';  % 相对于R_{i-1}坐标系，计算从R_{i-1}到第一个插值矩阵R_{i}的变换矩阵
    % end

    axang = rotm2axang(R(:, :, i));                %通过相邻两旋转矩阵间的轴角参数，相对于上一坐标系   
    axes(i, :) = axang(1:3);                %轴-角表示中的轴矢量
    angs(i)  = axang(4);                  %轴-角表示中的角位移标量,弧度
    %angles(i)  = rad2deg(axang(4));      %将弧度转换成角度

    % Calculate samll angles
    D_ang(i,:) = angs(i) * axes(i, :);

    % Calculate samll displacements
    D_e = [D_x D_y D_z D_ang(i,1) D_ang(i,2) D_ang(i,3)]';

    % Calculate the generalized inverse of Jacobian matrix
    Num_JM_L = JM_L(q(1),q(2),q(3),q(4),q(5),q(6),q(7));
    inv_JM_L = Num_JM_L'*(inv(Num_JM_L*Num_JM_L'));

    D_q(:,i) = inv_JM_L * D_e;
    q_Inter(:,i) = q + D_q(:,i);

    q = q_Inter(:,i);


    %Num_JM_L = JM_L(q_Inter(1,1),q_Inter(2,1),q_Inter(3,1),q_Inter(4,1),q_Inter(5,1),q_Inter(6,1),q_Inter(7,1))


end
figure
plot(q_Inter(1,:))
figure
plot(q_Inter(2,:))
figure
plot(q_Inter(3,:))
figure
plot(q_Inter(4,:))
figure
plot(q_Inter(5,:))
figure
plot(q_Inter(6,:))
figure
plot(q_Inter(7,:))



