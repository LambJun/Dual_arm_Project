clear;

syms x_A y_A z_A x_B y_B z_B x_C y_C z_C x_D y_D z_D x_E y_E z_E x_F y_F z_F x_G y_G z_G l_FG l_EG l_AD l_BC l_AB l_EF
assume(x_A,'real');
assume(y_A,'real');
assume(z_A,'real');
assume(x_B,'real');
assume(y_B,'real');
assume(z_B,'real');
assume(x_C,'real');
assume(y_C,'real');
assume(z_C,'real');
assume(x_D,'real');
assume(y_D,'real');
assume(z_D,'real');
assume(x_E,'real');
assume(y_E,'real');
assume(z_E,'real');
assume(x_F,'real');
assume(y_F,'real');
assume(z_F,'real');
assume(x_G,'real');
assume(y_G,'real');
assume(z_G,'real');
assume(l_FG,'positive');
assume(l_EG,'positive');
assume(l_AD,'positive');
assume(l_BC,'positive');
assume(l_AB,'positive');
assume(l_EF,'positive');


%AD
eq1 = (-x_B-x_D)^2 + (y_A-y_D)^2 + (z_A-z_D)^2 - l_AD^2;
%BC
eq2 = (x_B-x_C)^2 + (y_B-y_C)^2 + (z_B-z_C)^2 - l_BC^2;
%AB
eq3 = (-x_B-x_B)^2+ (y_A-y_B)^2 + (z_A-z_B)^2 - l_AB^2;
%AB perpendicular to EG
eq4 = [-x_B-x_B, y_A-y_B, z_A-z_B]*[(-x_B+x_B)/2-x_G, (y_A+y_B)/2-y_G, (z_A+z_B)/2-z_G]';
%EG
eq5 = (-x_G)^2 + ((y_A+y_B)/2-y_G)^2 + ((z_A+z_B)/2-z_G)^2 - l_EG^2;

f = [eq1, eq2, eq3, eq4, eq5];
variable = [x_B, y_B, z_B, y_A, z_A];
J = jacobian(f, variable);


% A and B initial position
x_B = 0.021;
y_B = 0.14123;
z_B = 0.01;
y_A = 0.14123;
z_A = 0.01;
% inputs
l_AD_0 = 0.1413;
% l_AD = 0.155;
l_BC_0 = 0.1413;
% l_BC = 0.16;
% other fixed variables
l_AB = 0.042;
l_EG = 0.015;
x_G = 0;
y_G = 0.14123;
z_G = -0.005;
x_D = -0.022;
y_D = 0;
z_D = 0.0025;
x_C = 0.022;
y_C = 0;
z_C = 0.0025;

x_n_plus_1 = [x_B, y_B, z_B, y_A, z_A]';% initialize
% disp([0, x_n_plus_1'])

figure(1)
% clf;
% x_A = -x_B;
% plot3(x_C,y_C,z_C,'r*');hold on;
% plot3(x_B,y_B,z_B,'g*');
% plot3([x_C,x_B,x_A,x_D],[y_C,y_B,y_A,y_D],[z_C,z_B,z_A,z_D]);
% plot3([(x_A+x_B)/2, x_G],[(y_A+y_B)/2, y_G],[(z_A+z_B)/2, z_G]);
% axis equal;

theta = 0:0.2:(2*pi);
centerX=0;
centerY=0.14123;
centerZ=-0.005;
R=0.015;
circle_y = R*cos(theta-1)+centerY;
circle_z = R*sin(theta-1)+centerZ;
for step = 1:10000
    disp("input")
    l_AD = 0.013 * sin(7*0.02*step) + 0.014 + 0.1189 + 0.0078%l_AD_0 + 0.014 * sin(7*0.02*step);
    l_BC = 0.013 * sin(13*0.002*step) + 0.014 + 0.1189 + 0.0078%l_BC_0 + 0.014 * sin(13*0.002*step);
    la_AD_pos = l_AD-0.1189-0.0078
    la_BC_pos = l_BC-0.1189-0.0078
    for ii = 1:100
        x_n = [x_B, y_B, z_B, y_A, z_A]';
        x_n_plus_1 = x_n - eval(J)\eval(f)';
        x_B = (x_n_plus_1(1)>0)*x_n_plus_1(1);
        y_B = (x_n_plus_1(2)>0)*x_n_plus_1(2);
        z_B = (x_n_plus_1(3)>0)*x_n_plus_1(3);
        y_A = (x_n_plus_1(4)>0)*x_n_plus_1(4);
        z_A = (x_n_plus_1(5)>0)*x_n_plus_1(5);
        if sum(abs(eval(f)))<1e-10
            % disp([step,ii])
            break;
        end
        if ii>8
            pause(0.01)
            disp([step,ii])
        end
    end

    %calculate Joint Position 
    GE0 = [0, 1];
    x_E = 0;
    y_E = (y_A+y_B)/2;
    z_E = (z_A+z_B)/2;
    GE = [y_E-centerY, z_E-centerZ];
    if y_E>y_G 
        joint6sign = 1;
    else
        joint6sign = -1;
    end
    disp("output")
    joint6pos = joint6sign * acos(GE*GE0' / (norm(GE) * norm(GE0)))/pi*180% rotate about x +
    EB0 = [1, 0, 0];
    EB = [x_B-x_E, y_B-y_E, z_B-z_E];
    if y_B>y_E
        joint7sign = 1;
    else
        joint7sign = -1;
    end
    joint7pos = joint7sign * acos(EB*EB0' / (norm(EB) * norm(EB0)))/pi*180% rotate about vector GE


    
    clf;
    plot3(theta*0,circle_y,circle_z);hold on;
    x_A = -x_B;
    plot3(x_C,y_C,z_C,'r*');
    plot3(x_B,y_B,z_B,'g*');
    plot3([x_C,x_B,x_A,x_D],[y_C,y_B,y_A,y_D],[z_C,z_B,z_A,z_D]);
    plot3([(x_A+x_B)/2, x_G],[(y_A+y_B)/2, y_G],[(z_A+z_B)/2, z_G]);
    axis equal;
    xlim([-0.03,0.03]);
    ylim([-0.001,0.17]);
    zlim([-0.01,0.025]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view([0.04,0.04,0.02])
    drawnow;
    pause(0.01);
end







