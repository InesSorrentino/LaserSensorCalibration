function [B_T_EEi, B_R_EEi, Si_q, B_D_i, Si_D_i] = sensingDataCalculation(data, rc)

% z axis of the sensor frame expressed in sensor frame
Si1_z_i = [0  0 1];

% Compute average chord direction expressed in sensor frame
Si1_c_i = 1/size(data.Si1_a_ij,1) * sum(data.Si1_b_ij - data.Si1_a_ij);

% Compute angle between the x axis of the sensor frame and the plane of the
% disk expressed in sensor frame
% theta = pi/2 - acos(dot(Si1_z_i, Si1_c_i) / (norm(Si1_z_i) * norm(Si1_c_i)));
theta = asin(abs(dot(Si1_z_i, Si1_c_i)) / (norm(Si1_z_i) * norm(Si1_c_i)));

R_yi_theta_neg = [cos(-theta), 0, sin(-theta);
  0         , 1,    0       ;
  -sin(-theta), 0, cos(-theta)];


% Rotate the end-points of the chords around yi of -theta
Si1r_a_ij = R_yi_theta_neg * data.Si1_a_ij';
Si1r_b_ij = R_yi_theta_neg * data.Si1_b_ij';

% Fit ellipse
Si1r_endpoints_x_ij(:,1) = [Si1r_a_ij(1,:)'; Si1r_b_ij(1,:)'];
Si1r_endpoints_z_ij(:,1) = [Si1r_a_ij(3,:)'; Si1r_b_ij(3,:)'];
[parameters, coefficients] = fitEllipse(Si1r_endpoints_x_ij, Si1r_endpoints_z_ij);

% Plot fitted ellipse
Cx = parameters(1);
Cz = parameters(2);
Rx = parameters(3);
Rz = parameters(4);
theta = parameters(5);
angle = linspace(0,2*pi);
X_ell = Rx * cos(angle);
Z_ell = Rz * sin(angle);
nx = X_ell * cos(theta) - Z_ell * sin(theta) + Cx;
nz = X_ell * sin(theta) + Z_ell * cos(theta) + Cz;
% figure,
% plot(nx, nz, 'o');
% xlabel('x');
% ylabel('z');
% title('Fitted ellipse in sensor frame');
% hold on;
% for ii = 1 : size(Si1r_b_ij,2)
%   line([Si1r_a_ij(1,ii), Si1r_b_ij(1,ii)], [Si1r_a_ij(3,ii), Si1r_b_ij(3,ii)], 'Color', 'red');
% end

a = coefficients(1);
b = coefficients(2);
c = coefficients(3);
d = coefficients(4);
e = coefficients(5);
f = coefficients(6);

% Compute five points ei in rotated sensor frame
Si1r_z1_tang = ( -(2*b*d - 4*a*e) - sqrt((2*b*d-4*a*e)^2 - 4*(b^2-4*a*c)*(d^2-4*a*f)) ) / ( 2*(b^2 - 4*a*c) );
Si1r_z2_tang = ( -(2*b*d - 4*a*e) + sqrt((2*b*d-4*a*e)^2 - 4*(b^2-4*a*c)*(d^2-4*a*f)) ) / ( 2*(b^2 - 4*a*c) );

Si1r_x1_tang = (-b*Si1r_z1_tang - d) / (2*a);
Si1r_x2_tang = (-b*Si1r_z2_tang - d) / (2*a);


Si1r_e1 = [Si1r_x1_tang; 0; Si1r_z1_tang];
Si1r_e2 = [Si1r_x2_tang; 0; Si1r_z2_tang];

Si1r_e3 = [parameters(1); 0; parameters(2)];

% Find lines e1e2, a1b1, ambm in the form y = alpha*x + beta in the plane
% x-z of the rotated sensor frame
alpha_e1e2 = (Si1r_e2(3) - Si1r_e1(3)) / (Si1r_e2(1) - Si1r_e1(1));
beta_e1e2 = - Si1r_e1(1)*(Si1r_e2(3) - Si1r_e1(3)) / (Si1r_e2(1) - Si1r_e1(1)) + Si1r_e1(3);

alpha_a1b1 = (Si1r_b_ij(3,1) - Si1r_a_ij(3,1)) / (Si1r_b_ij(1,1) - Si1r_a_ij(1,1));
beta_a1b1 = - Si1r_a_ij(1,1)*(Si1r_b_ij(3,1) - Si1r_a_ij(3,1)) / (Si1r_b_ij(1,1) - Si1r_a_ij(1,1)) + Si1r_a_ij(3,1);

alpha_ambm = (Si1r_b_ij(3,end) - Si1r_a_ij(3,end)) / (Si1r_b_ij(1,end) - Si1r_a_ij(1,end));
beta_ambm = - Si1r_a_ij(1,end)*(Si1r_b_ij(3,end) - Si1r_a_ij(3,end)) / (Si1r_b_ij(1,end) - Si1r_a_ij(1,end)) + Si1r_a_ij(3,end);

Si1r_e4(1,1) = (beta_a1b1 - beta_e1e2) / (alpha_e1e2 - alpha_a1b1);
Si1r_e4(2,1) = 0;
Si1r_e4(3,1) = alpha_e1e2*(beta_a1b1 - beta_e1e2)/(alpha_e1e2 - alpha_a1b1) + beta_e1e2;

Si1r_e5(1,1) = (beta_ambm - beta_e1e2) / (alpha_e1e2 - alpha_ambm);
Si1r_e5(2,1) = 0;
Si1r_e5(3,1) = alpha_e1e2*(beta_ambm - beta_e1e2)/(alpha_e1e2 - alpha_ambm) + beta_e1e2;

% scatter([Si1r_e1(1), Si1r_e2(1), Si1r_e3(1), Si1r_e4(1), Si1r_e5(1)], [Si1r_e1(3), Si1r_e2(3), Si1r_e3(3), Si1r_e4(3), Si1r_e5(3)]);
% plot([Si1r_e1(1), Si1r_e2(1), Si1r_e3(1), Si1r_e4(1), Si1r_e5(1)], [Si1r_e1(3), Si1r_e2(3), Si1r_e3(3), Si1r_e4(3), Si1r_e5(3)]);

% Compute from eta1 to eta5
etai1 = dot((Si1r_e1 - Si1r_e5), (Si1r_e4 - Si1r_e5)) / norm(Si1r_e4 - Si1r_e5)^2;
etai2 = dot((Si1r_e2 - Si1r_e5), (Si1r_e4 - Si1r_e5)) / norm(Si1r_e4 - Si1r_e5)^2;
etai3 = dot((Si1r_e3 - Si1r_e5), (Si1r_e4 - Si1r_e5)) / norm(Si1r_e4 - Si1r_e5)^2;
etai4 = dot((Si1r_e4 - Si1r_e5), (Si1r_e4 - Si1r_e5)) / norm(Si1r_e4 - Si1r_e5)^2;
etai5 = dot((Si1r_e5 - Si1r_e5), (Si1r_e4 - Si1r_e5)) / norm(Si1r_e4 - Si1r_e5)^2;


% Points e1, ..., e5 in the sensor frame
Si1_e1 = R_yi_theta_neg^(-1) * Si1r_e1;
Si1_e2 = R_yi_theta_neg^(-1) * Si1r_e2;
Si1_e3 = R_yi_theta_neg^(-1) * Si1r_e3;
Si1_e4 = R_yi_theta_neg^(-1) * Si1r_e4;
Si1_e5 = R_yi_theta_neg^(-1) * Si1r_e5;

etai11 = dot((Si1_e1 - Si1_e5), (Si1_e4 - Si1_e5)) / norm(Si1_e4 - Si1_e5)^2;
etai22 = dot((Si1_e2 - Si1_e5), (Si1_e4 - Si1_e5)) / norm(Si1_e4 - Si1_e5)^2;
etai33 = dot((Si1_e3 - Si1_e5), (Si1_e4 - Si1_e5)) / norm(Si1_e4 - Si1_e5)^2;
etai44 = dot((Si1_e4 - Si1_e5), (Si1_e4 - Si1_e5)) / norm(Si1_e4 - Si1_e5)^2;
etai55 = dot((Si1_e5 - Si1_e5), (Si1_e4 - Si1_e5)) / norm(Si1_e4 - Si1_e5)^2;


% Compute end-effector positions in base frame, corresponding to
% e1,...,e5 points
B_T_EEi1 = etai1 * data.kukaPose(1,1:3) + (1 - etai1) * data.kukaPose(end,1:3);
B_T_EEi2 = etai2 * data.kukaPose(1,1:3) + (1 - etai2) * data.kukaPose(end,1:3);
B_T_EEi3 = etai3 * data.kukaPose(1,1:3) + (1 - etai3) * data.kukaPose(end,1:3);
B_T_EEi4 = etai4 * data.kukaPose(1,1:3) + (1 - etai4) * data.kukaPose(end,1:3);
B_T_EEi5 = etai5 * data.kukaPose(1,1:3) + (1 - etai5) * data.kukaPose(end,1:3);

B_T_EEi = B_T_EEi3;

% Orientation is fixed and it is loaded from the dataset
B_abc_EEi = data.kukaPose(1,4:6);

B_abc_EEi = B_abc_EEi * pi / 180;

eul_A = B_abc_EEi(1);
eul_B = B_abc_EEi(2);
eul_C = B_abc_EEi(3);

Rx = [1,       0,              0;
     0,       cos(eul_C),   -sin(eul_C);
     0,       sin(eul_C),   cos(eul_C)];
   
Ry = [cos(eul_B),    0,      sin(eul_B);
      0,               1,      0;
     -sin(eul_B),   0,      cos(eul_B)];
   
Rz = [cos(eul_A),    -sin(eul_A),      0;
      sin(eul_A),    cos(eul_A),       0;
      0,               0,              1];
    
B_R_EEi = Rz*Ry*Rx;

% Position of the disk center in sensor frame
Si_q = R_yi_theta_neg^(-1) * Si1r_e3;

% End-effector offset in base frame
B_D_i = B_T_EEi2 - B_T_EEi1;

% Sensor offset in sensor frame
% Angle between the disk plane and the sensor plane
num = (norm(B_T_EEi2 - B_T_EEi1))^2 - (norm(Si1r_e1 - Si1r_e2))^2 - 4*rc^2;
den = 4*rc*(Si1r_e1(3) - Si1r_e2(3));
alphai = acos(num / den);

Si1r_D_i(1,1) = Si1r_e1(1) - Si1r_e2(1);
Si1r_D_i(2,1) = 2*rc*sin(alphai);
Si1r_D_i(3,1) = Si1r_e1(3) - Si1r_e2(3) + 2*rc*cos(alphai);

Si_D_i = R_yi_theta_neg^(-1) * Si1r_D_i;

if size(B_T_EEi,1) < 3
  B_T_EEi = B_T_EEi';
end
if size(Si_q,1) < 3
  Si_q = Si_q';
end
if size(B_D_i) < 3
  B_D_i = B_D_i';
end
if size(Si_D_i,1) < 3
  Si_D_i = Si_D_i';
end

end

