function [B_T_EEi, B_R_EEi, Si_q, B_D_i, Si_D_i] = sensingDataCalculation(data)

% z axis of the sensor frame expressed in sensor frame
Si1_z_i = [0  0 1];

% Compute average chord direction expressed in sensor frame
Si1_c_i = 1/size(data.Si1_a_ij,1) * sum(data.Si1_b_ij - data.Si1_a_ij);

% Compute angle between the x axis of the sensor frame and the plane of the
% disk expressed in sensor frame
theta = pi/2 - acos(abs(dot(Si1_z_i, Si1_c_i)) / (norm(Si1_z_i) * norm(Si1_c_i)));

% Rotate the end-points of the chords around yi of -theta
R_yi_theta_neg = [cos(-theta), 0, sin(-theta);
                   0         , 1,    0       ;
                 -sin(-theta), 0, cos(-theta)];


B_T_EEi = [];
B_R_EEi = [];
Si_q = [];
B_D_i = [];
Si_D_i = [];

end

