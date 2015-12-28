% Start simulating

% scanning parameters
delta_ti = 0.6;
slice_shifting_factor = 8;
actual_sampling_rate = delta_ti / slice_shifting_factor;

t = 0.04 : actual_sampling_rate : 7;
n_bolus = 7;

bolus_gap_factor = 1;
delta_bolus = bolus_gap_factor * delta_ti;


% parameters to be simulated
% reference values (to be searched)
cbf_0 = 60;
tau_0 = 0.6;
arrival_time = 0.7;
t1_tissue = 1.3;
t1_blood  = 1.6;

% simulate reference time series
asl_signal_ref = calculate_delta_M_tissue(cbf_0, tau_0, t, t1_tissue, t1_blood, arrival_time, n_bolus, delta_bolus);

% Grids
vector_cbf = 40 : 0.1 : 80;
vector_tau = 0.4 : 0.001 : 0.8;



cbf_size = length(vector_cbf);
tau_size = length(vector_tau);
ti_size  = length(t);

asl_signal_3D_matrix = zeros(cbf_size, tau_size, ti_size);
rmse_2D_matrix = zeros(cbf_size, tau_size);

for i = 1 : tau_size

	for j = 1 : cbf_size

		% generate time series here
		current_asl_signal = calculate_delta_M_tissue(vector_cbf(j), vector_tau(i), t, t1_tissue, t1_blood, arrival_time, n_bolus, delta_bolus);

		% assign this to signal matrix
		asl_signal_3D_matrix(j, i, :) = current_asl_signal(:);

		% calculate RMSE
		current_rmse = sqrt( mean((asl_signal_ref(:) - current_asl_signal(:)) .^ 2) );

		% assign to RMSE matrix
		rmse_2D_matrix(j, i) = current_rmse;

	end

end

another_siganl = asl_signal_3D_matrix(1, 1, :);
another_siganl = another_siganl(:);

figure;

plot(asl_signal_ref(:), 'b');

hold on;

plot(another_siganl(:), 'r');

figure;

surf(vector_tau, vector_cbf, rmse_2D_matrix);

local_mins = imregionalmin(rmse_2D_matrix, 4);

figure;
plot(local_mins);

% Grey scale image
figure;
imshow(mat2gray(local_mins));

