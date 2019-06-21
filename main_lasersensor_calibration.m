%% Load datasets
kukaPose_filename = 'KukaPose.txt';
chords_filename = 'endPoints.txt';

data1.folder = 'datasets/data1/segmented/';
data2.folder = 'datasets/data2/segmented/';
data3.folder = 'datasets/data3/segmented/';
data4.folder = 'datasets/data4/segmented/';
data5.folder = 'datasets/data5/segmented/';
data6.folder = 'datasets/data6/segmented/';
data7.folder = 'datasets/data7/segmented/';

data5.kukaPose = load([data5.folder, kukaPose_filename]);
data5.chords = load([data5.folder, chords_filename]);
data6.kukaPose = load([data6.folder, kukaPose_filename]);
data6.chords = load([data6.folder, chords_filename]);
data7.kukaPose = load([data7.folder, kukaPose_filename]);
data7.chords = load([data7.folder, chords_filename]);

data5.Si1_a_ij = [data5.chords(1:2:end,1), zeros(size(data5.chords,1)/2,1), data5.chords(1:2:end,2)];
data5.Si1_b_ij = [data5.chords(2:2:end,1), zeros(size(data5.chords,1)/2,1), data5.chords(2:2:end,2)];
data6.Si1_a_ij = [data6.chords(1:2:end,1), zeros(size(data5.chords,1)/2,1), data6.chords(1:2:end,2)];
data6.Si1_b_ij = [data6.chords(2:2:end,1), zeros(size(data5.chords,1)/2,1), data6.chords(2:2:end,2)];
data7.Si1_a_ij = [data7.chords(1:2:end,1), zeros(size(data5.chords,1)/2,1), data7.chords(1:2:end,2)];
data7.Si1_b_ij = [data7.chords(2:2:end,1), zeros(size(data5.chords,1)/2,1), data7.chords(2:2:end,2)];

%% Compute average chord direction
[data5.B_T_EEi, data5.B_R_EEi, data5.Si_q, data5.B_D_i, data5.Si_D_i] = sensingDataCalculation(data5);
[data6.B_T_EEi, data6.B_R_EEi, data6.Si_q, data6.B_D_i, data6.Si_D_i] = sensingDataCalculation(data6);
[data7.B_T_EEi, data7.B_R_EEi, data7.Si_q, data7.B_D_i, data7.Si_D_i] = sensingDataCalculation(data7);