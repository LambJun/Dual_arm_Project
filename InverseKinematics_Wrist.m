clear;
angle1=deg2rad(0);
posG = [0, 0, 0.13, 1];
theta = 0:0.2:(2*pi);

centerX=0; %G在坐标系中的位置
centerY=0.14088;
centerZ=-0.005;

R=0.015; %EG的长度
circle_y = R*cos(theta-2)+centerY;%E的位置
circle_z = R*sin(theta-2)+centerZ;

x_D=-0.022;
y_D=0;
z_D=0.0025;
x_C=0.022;
y_C=0;
z_C=0.0025;

LA_body_fixed_length = 0.1184+0.0005;

figure(1);
for ii = 1:10000
%% input
joint6position = sin(3*ii*0.01);
joint7posiiton = 0.23*sin(17*ii*0.001);
%% transform to A, B points. T6-rotx+rotx and T7-rotz
T6 = [1, 0, 0, 0;
      0, cos(joint6position), -sin(joint6position), 0;
      0, sin(joint6position), cos(joint6position), 0;
      0, 0, 0, 1]*...
      [1,0,0,0;
      0,1,0,0;
      0,0,1,0.015;
      0,0,0,1]*...
      [1, 0, 0, 0;
      0, cos(angle1), -sin(angle1), 0;
      0, sin(angle1), cos(angle1), 0;
      0, 0, 0, 1];%先转关节角；再沿z平移；再绕x轴旋转固定角度(这里的固定角度为偏置)。
T7A = [cos(joint7posiiton), -sin(joint7posiiton), 0, 0;
      sin(joint7posiiton), cos(joint7posiiton), 0, 0;
      0, 0, 1, 0.0;
      0, 0, 0, 1]*...
      [1,0,0,-0.021;
      0,1,0,0.0;
      0,0,1,0;
      0,0,0,1];
T7B = [cos(joint7posiiton), -sin(joint7posiiton), 0, 0;
      sin(joint7posiiton), cos(joint7posiiton), 0, 0;
      0, 0, 1, 0.0;
      0, 0, 0, 1]*...
      [1,0,0,0.021;
      0,1,0,0.0;
      0,0,1,0;
      0,0,0,1];
%% 
posA = T6*T7A;
posB = T6*T7B;
%% vector A and B
vectorA = [posA(1,4); posA(2,4); posA(3,4); 1]+[centerX,centerY,centerZ,0]';
vectorB = [posB(1,4); posB(2,4); posB(3,4); 1]+[centerX,centerY,centerZ,0]';
vectorE = [(vectorA(1)+vectorB(1))/2;(vectorA(2)+vectorB(2))/2;(vectorA(3)+vectorB(3))/2];

%%
clf
plot3(vectorA(1),vectorA(2),vectorA(3),'r*');hold on;
plot3(vectorB(1),vectorB(2),vectorB(3),'g*');
plot3(vectorE(1),vectorE(2),vectorE(3),'b*');
% disp((vectorA(1)+vectorB(1))/2)
% plot3([vectorA(1),vectorB(1)],[vectorA(2),vectorB(2)],[vectorA(3),vectorB(3)]);
% plot3([vectorA(1),0,vectorB(1)],[vectorA(2),0,vectorB(2)],[vectorA(3),0,vectorB(3)]);
plot3([centerX,vectorA(1),vectorE(1),vectorB(1),centerX], ...
        [centerY,vectorA(2),vectorE(2),vectorB(2),centerY], ...
        [centerZ,vectorA(3),vectorE(3),vectorB(3),centerZ] );
plot3([x_D,vectorA(1),vectorB(1),x_C], ...
        [y_D,vectorA(2),vectorB(2),y_C], ...
        [z_D,vectorA(3),vectorB(3),z_C],'k');
lengthAD=sqrt((vectorA(1)-x_D)^2+(vectorA(2)-y_D)^2+(vectorA(3)-z_D)^2);
lengthBC=sqrt((vectorB(1)-x_C)^2+(vectorB(2)-y_C)^2+(vectorB(3)-z_C)^2);
extendLengthAD=lengthAD-LA_body_fixed_length;
extendLengthBC=lengthBC-LA_body_fixed_length;
vectorAD=[x_D-vectorA(1),y_D-vectorA(2),z_D-vectorA(3)];
vectorBC=[x_C-vectorB(1),y_C-vectorB(2),z_C-vectorB(3)];
plot3([vectorA(1),vectorA(1)+vectorAD(1)/lengthAD*extendLengthAD], ...
        [vectorA(2),vectorA(2)+vectorAD(2)/lengthAD*extendLengthAD], ...
        [vectorA(3),vectorA(3)+vectorAD(3)/lengthAD*extendLengthAD],'r');
plot3([vectorB(1),vectorB(1)+vectorBC(1)/lengthBC*extendLengthBC], ...
        [vectorB(2),vectorB(2)+vectorBC(2)/lengthBC*extendLengthBC], ...
        [vectorB(3),vectorB(3)+vectorBC(3)/lengthBC*extendLengthBC],'g');
disp([joint6position,joint7posiiton])
% disp([extendLengthAD,extendLengthBC])
% x_E = T6(1,4); y_E = T6(2,4); z_E = T6(3,4); 
plot3([vectorE(1),vectorE(1)+T6(1,3)*0.01],[vectorE(2),vectorE(2)+T6(2,3)*0.01],[vectorE(3),vectorE(3)+T6(3,3)*0.01]);

plot3(theta*0,circle_y,circle_z);

axis equal;
xlim([-0.03,0.03]);
ylim([-0.001,0.17]);
zlim([-0.01,0.025]);
xlabel('x');
ylabel('y');
zlabel('z');
view([0.04,0.04,0.02])
drawnow;
pause(0.01)
end


