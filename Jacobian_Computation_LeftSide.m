
% Standard DH model (SDH) is adopted for modeling
% This is for bi-manipulator Project
% Created on 11/11/2024 
clc;
clear;

syms theta_L1 theta_L2 theta_L3 theta_L4 theta_L5 theta_L6 theta_L7 real
syms theta_R1 theta_R2 theta_R3 theta_R4 theta_R5 theta_R6 theta_R7 real

theta_L = sym('theta_L',[7,1], 'real');
theta_R = sym('theta_R',[7,1], 'real');

A_L = sym('AL',[4,4], 'real');
A_R = sym('AR',[4,4], 'real');

A_0E_L = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A_0E_R = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];

Z  = sym('Z',[3,6], 'real');
P  = sym('P',[3,7], 'real');
PT = sym('PT',[4,1], 'real');
Jp = sym('Jp',[3,7], 'real');
Jo = sym('Jo',[3,7], 'real');
JL_SG = sym('JL_SG',[6,7], 'real');

% Z_0 = [0; 0; 1];
% P_0 = [0; 0.121; 0.071];
% EP_0 = [0; 0; 0; 1];

Z_0 = [0; 0; 1];
P_0 = [0; 0; 0];
EP_0 = [0; 0; 0; 1];

TheV = 1.0000e-10; % Threshold for simplifying symbolic expressions


%% DH parameters derived by SDH model

% DH parameters of left arm
%        theta（angle）      d              a              alpha（offset）         offset(theta,original position)
left(1,:) = [0.0           0.084           0.0                 pi/2                 pi/2 ]; 
left(2,:) = [0.0           0.0             0.0                 pi/2                 pi/2 ]; 
left(3,:) = [0.0           0.27            0.0                 pi/2                 pi/2 ]; 
left(4,:) = [0.0           0.0             0.0                 pi/2                 pi   ]; 
left(5,:) = [0.0           0.275           0.005               pi/2                 -pi/2]; 
left(6,:) = [0.0           0.0             0.0                 pi/2                 pi/2 ]; 
left(7,:) = [0.0           0.0             0.023               0                    0    ]; 

% DH parameters of right arm
%         theta（angle）     d              a              alpha（offset）         offset(theta,original postion)
right(1,:) = [0.0           0.084          0.0                 pi/2                 pi/2 ]; 
right(2,:) = [0.0           0.0            0.0                 pi/2                 -pi/2]; 
right(3,:) = [0.0           0.27           0.0                 pi/2                 -pi/2]; 
right(4,:) = [0.0           0.0            0.0                 pi/2                 pi   ]; 
right(5,:) = [0.0           0.275          0.005               -pi/2                -pi/2]; 
right(6,:) = [0.0           0.0            0.0                 pi/2                 pi/2 ]; 
right(7,:) = [0.0           0.0            0.023               0                    pi   ]; 

%% Calculations for the left arm

%% Homogeneous transformation matrix of the left arm

% Homogeneous transformation from the Zero coordinate to the End coordinate 
for i=1:7

    % Transformation matrix from Frame i to Frame i-1
    A_L(:,:) = [cos(theta_L(i)+left(i,5)),   -sin(theta_L(i)+left(i,5))*cos(left(i,4)),    sin(theta_L(i)+left(i,5))*sin(left(i,4)),   left(i,3)*cos(theta_L(i)+left(i,5));
                sin(theta_L(i)+left(i,5)),    cos(theta_L(i)+left(i,5))*cos(left(i,4)),   -cos(theta_L(i)+left(i,5))*sin(left(i,4)),   left(i,3)*sin(theta_L(i)+left(i,5));
                0,                            sin(left(i,4)),                              cos(left(i,4)),                             left(i,2);
                0,                            0,                                           0,                                          1                                   ];

    % Transformation matrix recursion
    A_0E_L = A_0E_L*A_L; 
    % Calculate Z_{i}
    if i<7
      Z(1:3,i)= A_0E_L(1:3,1:3)*Z_0;
    end
    % Calculate P_{i}
      PT= A_0E_L*EP_0;
      P(1:3,i)= PT(1:3,1);

end

% Homogeneous transformation from the Global coordinate to the Zero coordinate
A_G0_L = [1, 0, 0, 0; 0, 0, 1, 0.121; 0, -1, 0, 0.071; 0, 0, 0, 1];
R_G0_L = [1, 0, 0; 0, 0, 1; 0, -1, 0];

% Homogeneous transformation from the Global coordinate to the End coordinate
A_GE_L = A_G0_L * A_0E_L;

%% Calculate the Jacobian Matrix of the lef arm
% Calculate the first colum of the Jacobian Matrix

     P_e = P(1:3,7)-P_0;
     Jp(1:3,1) = [Z_0(2)*P_e(3)-Z_0(3)*P_e(2); -Z_0(1)*P_e(3)+Z_0(3)*P_e(1); Z_0(1)*P_e(2)-Z_0(2)*P_e(1)];
     Jo(1:3,1) = Z_0;
% Calculate the last six colums of the Jacobian Matrix
 for i=1:6

     P_e = P(1:3,7)-P(1:3,i);
     Jp(1:3,i+1) = [Z(2,i)*P_e(3)-Z(3,i)*P_e(2); -Z(1,i)*P_e(3)+Z(3,i)*P_e(1); Z(1,i)*P_e(2)-Z(2,i)*P_e(1)];
     Jo(1:3,i+1) = Z(1:3,i);        

end

% Resulting in Jacobian Matrix of the left arm within the base coordinate
    JL = [Jp;Jo];

% Jacobian Matrix within the global coordinate
    digits(10) %             % Save ten effective digits
    R_GO = [R_G0_L, zeros(3,3); zeros(3,3), R_G0_L];
    JL_G = R_GO* JL;         % Jacobian Matrix of the left side in the global coordinate
    JL_G_1 = vpa(JL_G,10);   % Simplify the expression of coefficients and the number of digits
    JL_G_2 = expand(JL_G_1); % Simplify the expression of JL_G_1 in the form of multiplication of coefficients and symbolic functions

    digits(10)
    J_0= real(vpa(subs(JL_G,[theta_L1,theta_L2,theta_L3,theta_L4,theta_L5,theta_L6,theta_L7],[0,0,pi/2,pi/3,0,pi/9,0])));
    J_1= real(vpa(subs(JL_G_1,[theta_L1,theta_L2,theta_L3,theta_L4,theta_L5,theta_L6,theta_L7],[0,0,pi/2,pi/3,0,pi/9,0])));
    J_2= real(vpa(subs(JL_G_2,[theta_L1,theta_L2,theta_L3,theta_L4,theta_L5,theta_L6,theta_L7],[0,0,pi/2,pi/3,0,pi/9,0])));
    

% Separate the coefficients and symbolic functions of all elements in the first row within the Jacobian Matrix    
    [coeff_L11, terms_L11] = coeffs(JL_G_2(1,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L12, terms_L12] = coeffs(JL_G_2(1,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L13, terms_L13] = coeffs(JL_G_2(1,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L14, terms_L14] = coeffs(JL_G_2(1,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L15, terms_L15] = coeffs(JL_G_2(1,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L16, terms_L16] = coeffs(JL_G_2(1,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L17, terms_L17] = coeffs(JL_G_2(1,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);

% Separate the coefficients and symbolic functions of all elements in the second row within the Jacobian Matrix    
    [coeff_L21, terms_L21] = coeffs(JL_G_2(2,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L22, terms_L22] = coeffs(JL_G_2(2,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L23, terms_L23] = coeffs(JL_G_2(2,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L24, terms_L24] = coeffs(JL_G_2(2,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L25, terms_L25] = coeffs(JL_G_2(2,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L26, terms_L26] = coeffs(JL_G_2(2,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L27, terms_L27] = coeffs(JL_G_2(2,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);

% Separate the coefficients and symbolic functions of all elements in the third row within the Jacobian Matrix    
    [coeff_L31, terms_L31] = coeffs(JL_G_2(3,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L32, terms_L32] = coeffs(JL_G_2(3,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L33, terms_L33] = coeffs(JL_G_2(3,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L34, terms_L34] = coeffs(JL_G_2(3,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L35, terms_L35] = coeffs(JL_G_2(3,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L36, terms_L36] = coeffs(JL_G_2(3,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L37, terms_L37] = coeffs(JL_G_2(3,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);

% Separate the coefficients and symbolic functions of all elements in the fourth row within the Jacobian Matrix    
    [coeff_L41, terms_L41] = coeffs(JL_G_2(4,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L42, terms_L42] = coeffs(JL_G_2(4,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L43, terms_L43] = coeffs(JL_G_2(4,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L44, terms_L44] = coeffs(JL_G_2(4,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L45, terms_L45] = coeffs(JL_G_2(4,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L46, terms_L46] = coeffs(JL_G_2(4,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L47, terms_L47] = coeffs(JL_G_2(4,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);

% Separate the coefficients and symbolic functions of all elements in the fifth row within the Jacobian Matrix    
    [coeff_L51, terms_L51] = coeffs(JL_G_2(5,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L52, terms_L52] = coeffs(JL_G_2(5,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L53, terms_L53] = coeffs(JL_G_2(5,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L54, terms_L54] = coeffs(JL_G_2(5,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L55, terms_L55] = coeffs(JL_G_2(5,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L56, terms_L56] = coeffs(JL_G_2(5,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L57, terms_L57] = coeffs(JL_G_2(5,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);

% Separate the coefficients and symbolic functions of all elements in the sixth row within the Jacobian Matrix    
    [coeff_L61, terms_L61] = coeffs(JL_G_2(6,1),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L62, terms_L62] = coeffs(JL_G_2(6,2),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L63, terms_L63] = coeffs(JL_G_2(6,3),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L64, terms_L64] = coeffs(JL_G_2(6,4),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L65, terms_L65] = coeffs(JL_G_2(6,5),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L66, terms_L66] = coeffs(JL_G_2(6,6),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
    [coeff_L67, terms_L67] = coeffs(JL_G_2(6,7),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);


% Simplify the symbolic expression of each element by setting the
% coefficients with absolute values less than 1e-10 to zero

    % Get the Maximum dimension length of each coefficient vector
    n_L11=length(coeff_L11); n_L12=length(coeff_L12); n_L13=length(coeff_L13); n_L14=length(coeff_L14); n_L15=length(coeff_L15); n_L16=length(coeff_L16); n_L17=length(coeff_L17); 
    n_L21=length(coeff_L21); n_L22=length(coeff_L22); n_L23=length(coeff_L23); n_L24=length(coeff_L24); n_L25=length(coeff_L25); n_L26=length(coeff_L26); n_L27=length(coeff_L27); 
    n_L31=length(coeff_L31); n_L32=length(coeff_L32); n_L33=length(coeff_L33); n_L34=length(coeff_L34); n_L35=length(coeff_L35); n_L36=length(coeff_L36); n_L37=length(coeff_L37); 
    n_L41=length(coeff_L41); n_L42=length(coeff_L42); n_L43=length(coeff_L43); n_L44=length(coeff_L44); n_L45=length(coeff_L45); n_L46=length(coeff_L46); n_L47=length(coeff_L47); 
    n_L51=length(coeff_L51); n_L52=length(coeff_L52); n_L53=length(coeff_L53); n_L54=length(coeff_L54); n_L55=length(coeff_L55); n_L56=length(coeff_L56); n_L57=length(coeff_L57); 
    n_L61=length(coeff_L61); n_L62=length(coeff_L62); n_L63=length(coeff_L63); n_L64=length(coeff_L64); n_L65=length(coeff_L65); n_L66=length(coeff_L66); n_L67=length(coeff_L67); 

    n_L = [n_L11, n_L12, n_L13, n_L14, n_L15, n_L16, n_L17;
           n_L21, n_L22, n_L23, n_L24, n_L25, n_L26, n_L27;
           n_L31, n_L32, n_L33, n_L34, n_L35, n_L36, n_L37;
           n_L41, n_L42, n_L43, n_L44, n_L45, n_L46, n_L47;
           n_L51, n_L52, n_L53, n_L54, n_L55, n_L56, n_L57;
           n_L61, n_L62, n_L63, n_L64, n_L65, n_L66, n_L67 ];

    % Define an internediate variable that can be used to compare coefficient vectors before and after simplification
    coeff_L11_1=coeff_L11;  coeff_L12_1=coeff_L12; coeff_L13_1=coeff_L13; coeff_L14_1=coeff_L14; coeff_L15_1=coeff_L15; coeff_L16_1=coeff_L16; coeff_L17_1=coeff_L17;
    coeff_L21_1=coeff_L21;  coeff_L22_1=coeff_L22; coeff_L23_1=coeff_L23; coeff_L24_1=coeff_L24; coeff_L25_1=coeff_L25; coeff_L26_1=coeff_L26; coeff_L27_1=coeff_L27;
    coeff_L31_1=coeff_L31;  coeff_L32_1=coeff_L32; coeff_L33_1=coeff_L33; coeff_L34_1=coeff_L34; coeff_L35_1=coeff_L35; coeff_L36_1=coeff_L36; coeff_L37_1=coeff_L37;
    coeff_L41_1=coeff_L41;  coeff_L42_1=coeff_L42; coeff_L43_1=coeff_L43; coeff_L44_1=coeff_L44; coeff_L45_1=coeff_L45; coeff_L46_1=coeff_L46; coeff_L47_1=coeff_L47;
    coeff_L51_1=coeff_L51;  coeff_L52_1=coeff_L52; coeff_L53_1=coeff_L53; coeff_L54_1=coeff_L54; coeff_L55_1=coeff_L55; coeff_L56_1=coeff_L56; coeff_L57_1=coeff_L57;
    coeff_L61_1=coeff_L61;  coeff_L62_1=coeff_L62; coeff_L63_1=coeff_L63; coeff_L64_1=coeff_L64; coeff_L65_1=coeff_L65; coeff_L66_1=coeff_L66; coeff_L67_1=coeff_L67;

    % Processing JL_G_2(1,1) within Jacobian matrix  
    for i=1:n_L(1,1)

        if abs(coeff_L11_1(i)) < TheV % if the value is less than the threshold value
            coeff_L11_1(i)=0;
        end

    end
    % Processing JL_G_2(1,2) within Jacobian matrix
    for i=1:n_L(1,2)

        if abs(coeff_L12_1(i)) < TheV
            coeff_L12_1(i)=0;
        end

    end
    % Processing JL_G_2(1,3) within Jacobian matrix
    for i=1:n_L(1,3)

        if abs(coeff_L13_1(i)) < TheV 
            coeff_L13_1(i)=0;
        end

    end
    % Processing JL_G_2(1,4) within Jacobian matrix
    for i=1:n_L(1,4)

        if abs(coeff_L14_1(i)) < TheV 
            coeff_L14_1(i)=0;
        end

    end
    % Processing JL_G_2(1,5) within Jacobian matrix
    for i=1:n_L(1,5)

        if abs(coeff_L15_1(i)) < TheV 
            coeff_L15_1(i)=0;
        end

    end
    % Processing JL_G_2(1,6) within Jacobian matrix
    for i=1:n_L(1,6)

        if abs(coeff_L16_1(i)) < TheV 
            coeff_L16_1(i)=0;
        end

    end
    % Processing JL_G_2(1,7) within Jacobian matrix
    for i=1:n_L(1,7)

        if abs(coeff_L17_1(i)) < TheV 
            coeff_L17_1(i)=0;
        end

    end

    % Processing JL_G_2(2,1) within Jacobian matrix
    for i=1:n_L(2,1)

        if abs(coeff_L21_1(i)) < TheV
            coeff_L21_1(i)=0;
        end

    end 
    % Processing JL_G_2(2,2) within Jacobian matrix
    for i=1:n_L(2,2)

        if abs(coeff_L22_1(i)) < TheV 
            coeff_L22_1(i)=0;
        end

    end
    % Processing JL_G_2(2,3) within Jacobian matrix
    for i=1:n_L(2,3)

        if abs(coeff_L23_1(i)) < TheV 
            coeff_L23_1(i)=0;
        end

    end
    % Processing JL_G_2(2,4) within Jacobian matrix
    for i=1:n_L(2,4)

        if abs(coeff_L24_1(i)) < TheV 
            coeff_L24_1(i)=0;
        end

    end
    % Processing JL_G_2(2,5) within Jacobian matrix
    for i=1:n_L(2,5)

        if abs(coeff_L25_1(i)) < TheV 
            coeff_L25_1(i)=0;
        end

    end
    % Processing JL_G_2(2,6) within Jacobian matrix
    for i=1:n_L(2,6)

        if abs(coeff_L26_1(i)) < TheV 
            coeff_L26_1(i)=0;
        end

    end
    % Processing JL_G_2(2,7) within Jacobian matrix
    for i=1:n_L(2,7)

        if abs(coeff_L27_1(i)) < TheV 
            coeff_L27_1(i)=0;
        end

    end

    % Processing JL_G_2(3,1) within Jacobian matrix
    for i=1:n_L(3,1)

        if abs(coeff_L31_1(i)) < TheV
            coeff_L31_1(i)=0;
        end

    end
    % Processing JL_G_2(3,2) within Jacobian matrix
    for i=1:n_L(3,2)

        if abs(coeff_L32_1(i)) < TheV
            coeff_L32_1(i)=0;
        end

    end
    % Processing JL_G_2(3,3) within Jacobian matrix
    for i=1:n_L(3,3)

        if abs(coeff_L33_1(i)) < TheV 
            coeff_L33_1(i)=0;
        end

    end
    % Processing JL_G_2(3,4) within Jacobian matrix
    for i=1:n_L(3,4)

        if abs(coeff_L34_1(i)) < TheV 
            coeff_L34_1(i)=0;
        end

    end
    % Processing JL_G_2(3,5) within Jacobian matrix
    for i=1:n_L(3,5)

        if abs(coeff_L35_1(i)) < TheV 
            coeff_L35_1(i)=0;
        end

    end
    % Processing JL_G_2(3,6) within Jacobian matrix
    for i=1:n_L(3,6)

        if abs(coeff_L36_1(i)) < TheV 
            coeff_L36_1(i)=0;
        end

    end
    % Processing JL_G_2(3,7) within Jacobian matrix
    for i=1:n_L(3,7)

        if abs(coeff_L37_1(i)) < TheV 
            coeff_L37_1(i)=0;
        end

    end

    % Processing JL_G_2(4,1) within Jacobian matrix
    for i=1:n_L(4,1)

        if abs(coeff_L41_1(i)) < TheV
            coeff_L41_1(i)=0;
        end

    end
    % Processing JL_G_2(4,2) within Jacobian matrix
    for i=1:n_L(4,2)

        if abs(coeff_L42_1(i)) < TheV 
            coeff_L42_1(i)=0;
        end

    end
    % Processing JL_G_2(4,3) within Jacobian matrix
    for i=1:n_L(4,3)

        if abs(coeff_L43_1(i)) < TheV 
            coeff_L43_1(i)=0;
        end

    end
    % Processing JL_G_2(4,4) within Jacobian matrix
    for i=1:n_L(4,4)

        if abs(coeff_L44_1(i)) < TheV 
            coeff_L44_1(i)=0;
        end

    end
    % Processing JL_G_2(4,5) within Jacobian matrix
    for i=1:n_L(4,5)

        if abs(coeff_L45_1(i)) < TheV 
            coeff_L45_1(i)=0;
        end

    end
    % Processing JL_G_2(4,6) within Jacobian matrix
    for i=1:n_L(4,6)

        if abs(coeff_L46_1(i)) < TheV 
            coeff_L46_1(i)=0;
        end

    end
    % Processing JL_G_2(4,7) within Jacobian matrix
    for i=1:n_L(4,7)

        if abs(coeff_L47_1(i)) < TheV 
            coeff_L47_1(i)=0;
        end

    end

    % Processing JL_G_2(5,1) within Jacobian matrix
    for i=1:n_L(5,1)

        if abs(coeff_L51_1(i)) < TheV
            coeff_L51_1(i)=0;
        end

    end
    % Processing JL_G_2(5,2) within Jacobian matrix
    for i=1:n_L(5,2)

        if abs(coeff_L52_1(i)) < TheV
            coeff_L52_1(i)=0;
        end

    end
    % Processing JL_G_2(5,3) within Jacobian matrix
    for i=1:n_L(5,3)

        if abs(coeff_L53_1(i)) < TheV 
            coeff_L53_1(i)=0;
        end

    end
    % Processing JL_G_2(5,4) within Jacobian matrix
    for i=1:n_L(5,4)

        if abs(coeff_L54_1(i)) < TheV 
            coeff_L54_1(i)=0;
        end

    end
    % Processing JL_G_2(5,5) within Jacobian matrix
    for i=1:n_L(5,5)

        if abs(coeff_L55_1(i)) < TheV 
            coeff_L55_1(i)=0;
        end

    end
    % Processing JL_G_2(5,6) within Jacobian matrix
    for i=1:n_L(5,6)

        if abs(coeff_L56_1(i)) < TheV 
            coeff_L56_1(i)=0;
        end

    end
    % Processing JL_G_2(5,7) within Jacobian matrix
    for i=1:n_L(5,7)

        if abs(coeff_L57_1(i)) < TheV 
            coeff_L57_1(i)=0;
        end

    end

    % Processing JL_G_2(6,1) within Jacobian matrix
    for i=1:n_L(6,1)

        if abs(coeff_L61_1(i)) < TheV
            coeff_L61_1(i)=0;
        end

    end
    % Processing JL_G_2(6,2) within Jacobian matrix
    for i=1:n_L(6,2)

        if abs(coeff_L62_1(i)) < TheV
            coeff_L62_1(i)=0;
        end

    end
    % Processing JL_G_2(6,3) within Jacobian matrix
    for i=1:n_L(6,3)

        if abs(coeff_L63_1(i)) < TheV 
            coeff_L63_1(i)=0;
        end

    end
    % Processing JL_G_2(6,4) within Jacobian matrix
    for i=1:n_L(6,4)

        if abs(coeff_L64_1(i)) < TheV 
            coeff_L64_1(i)=0;
        end

    end
    % Processing JL_G_2(6,5) within Jacobian matrix
    for i=1:n_L(6,5)

        if abs(coeff_L65_1(i)) < TheV 
            coeff_L65_1(i)=0;
        end

    end
    % Processing JL_G_2(6,6) within Jacobian matrix
    for i=1:n_L(6,6)

        if abs(coeff_L66_1(i)) < TheV 
            coeff_L66_1(i)=0;
        end

    end
    % Processing JL_G_2(6,7) within Jacobian matrix
    for i=1:n_L(6,7)

        if abs(coeff_L67_1(i)) < TheV 
            coeff_L67_1(i)=0;
        end

    end

    % Simplified symbolic expression of Jacobian matrix within the global coordinate
    JL_SG(1,1)=coeff_L11_1*terms_L11'; JL_SG(1,2)=coeff_L12_1*terms_L12'; JL_SG(1,3)=coeff_L13_1*terms_L13'; JL_SG(1,4)=coeff_L14_1*terms_L14'; JL_SG(1,5)=coeff_L15_1*terms_L15'; JL_SG(1,6)=coeff_L16_1*terms_L16'; JL_SG(1,7)=coeff_L17_1*terms_L17'; 
    JL_SG(2,1)=coeff_L21_1*terms_L21'; JL_SG(2,2)=coeff_L22_1*terms_L22'; JL_SG(2,3)=coeff_L23_1*terms_L23'; JL_SG(2,4)=coeff_L24_1*terms_L24'; JL_SG(2,5)=coeff_L25_1*terms_L25'; JL_SG(2,6)=coeff_L26_1*terms_L26'; JL_SG(2,7)=coeff_L27_1*terms_L27';
    JL_SG(3,1)=coeff_L31_1*terms_L31'; JL_SG(3,2)=coeff_L32_1*terms_L32'; JL_SG(3,3)=coeff_L33_1*terms_L33'; JL_SG(3,4)=coeff_L34_1*terms_L34'; JL_SG(3,5)=coeff_L35_1*terms_L35'; JL_SG(3,6)=coeff_L36_1*terms_L36'; JL_SG(3,7)=coeff_L37_1*terms_L37'; 
    JL_SG(4,1)=coeff_L41_1*terms_L41'; JL_SG(4,2)=coeff_L42_1*terms_L42'; JL_SG(4,3)=coeff_L43_1*terms_L43'; JL_SG(4,4)=coeff_L44_1*terms_L44'; JL_SG(4,5)=coeff_L45_1*terms_L45'; JL_SG(4,6)=coeff_L46_1*terms_L46'; JL_SG(4,7)=coeff_L47_1*terms_L47';
    JL_SG(5,1)=coeff_L51_1*terms_L51'; JL_SG(5,2)=coeff_L52_1*terms_L52'; JL_SG(5,3)=coeff_L53_1*terms_L53'; JL_SG(5,4)=coeff_L54_1*terms_L54'; JL_SG(5,5)=coeff_L55_1*terms_L55'; JL_SG(5,6)=coeff_L56_1*terms_L56'; JL_SG(5,7)=coeff_L57_1*terms_L57'; 
    JL_SG(6,1)=coeff_L61_1*terms_L61'; JL_SG(6,2)=coeff_L62_1*terms_L62'; JL_SG(6,3)=coeff_L63_1*terms_L63'; JL_SG(6,4)=coeff_L64_1*terms_L64'; JL_SG(6,5)=coeff_L65_1*terms_L65'; JL_SG(6,6)=coeff_L66_1*terms_L66'; JL_SG(6,7)=coeff_L67_1*terms_L67';

    % Get the number of contributor of each element within simplified Jacobian matrix
    n_SL = zeros(6,7);
    for i=1:6
       for j=1:7

         [coeff_SL, terms_SL] = coeffs(JL_SG(i,j),[cos(theta_L1),cos(theta_L2),cos(theta_L3),cos(theta_L4),cos(theta_L5),cos(theta_L6),cos(theta_L7),sin(theta_L1),sin(theta_L2),sin(theta_L3),sin(theta_L4),sin(theta_L5),sin(theta_L6),sin(theta_L7)]);
         n_SL(i,j) = length(coeff_SL); % Get the number of contributor of each element within simplified Jacobian matrix

       end
    end

    J_3= real(vpa(subs(JL_SG,[theta_L1,theta_L2,theta_L3,theta_L4,theta_L5,theta_L6,theta_L7],[0,0,pi/2,pi/3,0,pi/9,0])));

    % 
    matlabFunction(JL_SG, 'File', 'JM_L');









