function [B_T_EEi, B_R_EEi, Si_q, B_D_i, Si_D_i] = sensingDataCalculation(data)

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
figure,
plot(nx, nz, 'o');
xlabel('x');
ylabel('z');
title('Fitted ellipse in sensor frame');
hold on;
for ii = 1 : size(Si1r_b_ij,2)
  line([Si1r_a_ij(1,ii), Si1r_b_ij(1,ii)], [Si1r_a_ij(3,ii), Si1r_b_ij(3,ii)], 'Color', 'red');
end

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


% Inverti punti
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

scatter(Si1r_e4(1),Si1r_e4(3));
scatter(Si1r_e5(1),Si1r_e5(3));
plot([Si1r_e1(1), Si1r_e2(1)], [Si1r_e1(3), Si1r_e2(3)]);




B_T_EEi = [];
B_R_EEi = [];
Si_q = [];
B_D_i = [];
Si_D_i = [];

end

