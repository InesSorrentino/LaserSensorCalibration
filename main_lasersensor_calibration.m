%% Load datasets
kukaPose_filename = 'KukaPose.txt';
chords_filename = 'endPoints.txt';

data(1).folder = 'datasets/data1/segmented/';
data(2).folder = 'datasets/data2/segmented/';
data(3).folder = 'datasets/data3/segmented/';
data(4).folder = 'datasets/data4/segmented/';
data(5).folder = 'datasets/data5/segmented/';
data(6).folder = 'datasets/data6/segmented/';
data(7).folder = 'datasets/data7/segmented/';

for i = 1 : 7
  if i ~= 2
    data(i).kukaPose = load([data(i).folder, kukaPose_filename]);
    data(i).chords = load([data(i).folder, chords_filename]);
  end
end

for i = 1 : 7
  if i ~= 2
    data(i).Si1_a_ij = [data(i).chords(1:2:end,1), zeros(size(data(i).chords,1)/2,1), data(i).chords(1:2:end,2)];
    data(i).Si1_b_ij = [data(i).chords(2:2:end,1), zeros(size(data(i).chords,1)/2,1), data(i).chords(2:2:end,2)];
  end
end


%% Compute average chord direction
rc = 80;

for i = 1 : 7
  if i ~= 2
    [data(i).B_T_EEi(:,1), data(i).B_R_EEi, data(i).Si_q(:,1), data(i).B_D_i(:,1), data(i).Si_D_i(:,1)] = sensingDataCalculation(data(i), rc);
  end
end


indeces_to_use = [1 5 7];
X = [];
Y = [];

for i = 1 : length(indeces_to_use)
  X = [X, data(indeces_to_use(i)).Si_D_i];
  Y = [Y, data(indeces_to_use(i)).B_R_EEi^(-1)*data(indeces_to_use(i)).B_D_i];
end

H = X * Y';

[U,S,V] = svd(H);

EE_R_S = V * U';


M = [];
N = [];

for i = 1 : (length(indeces_to_use)-1)
  M = [M; data(indeces_to_use(i)).B_R_EEi - data(indeces_to_use(i+1)).B_R_EEi];
  N = [N; data(indeces_to_use(i+1)).B_R_EEi*EE_R_S*data(indeces_to_use(i+1)).Si_q - ...
          data(indeces_to_use(i)).B_R_EEi*EE_R_S*data(indeces_to_use(i)).Si_q + ...
          data(indeces_to_use(i+1)).B_T_EEi - ...
          data(indeces_to_use(i)).B_T_EEi];
end


EE_T_S = (M' * M)^(-1) * M' * N;
 
