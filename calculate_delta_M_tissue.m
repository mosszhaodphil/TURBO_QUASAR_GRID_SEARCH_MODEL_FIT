% This code is the implementation the following papers
% MA Chappell (2012) doi: 10.1002/mrm.24372 (MACQ)
% MA Chappell (2012) doi: 10.1002/mrm.24260 (MACD)
% ET Petersen (2006) doi: 10.1002/mrm.20784 (ETP)
% RB Buxton (1998) doi: 10.1002/mrm.1910400308 (RBB)
% L Ostergaard (1996) doi: 10.1002/mrm.1910360510 (LO)
% J Hrabe (2003) doi:10.1016/j.jmr.2003.11.002 (JH)

% This function calculates ASL signal deltaM of tissue using Buxton's model (RBB)
% The same method is also used in equation [5] of (MACQ)
% delta_M_tissue = 2 * alpha * M0a * f * (c(t) * r(t) * m(t))
% c(t) = exp(-1 / T1a) * a(t)
% a(t) depends on dispersion

% Here we adopt JH's analytical solution to calculate tissue magnetization

function delta_M_tissue = calculate_delta_M_tissue(cbf, tau, t, t1_tissue, t1_blood, arrival_time, n_bolus, delta_bolus)

	delta_M_tissue = zeros(length(t), 1); % ASL signal of tissue

	bolus_arrived = 0;

	% Check dispersion
	dispersion_type = 1;

	% No dispersion
	if(dispersion_type == 1)
		while(bolus_arrived < n_bolus)

			%bolus_time_passed = bolus_arrived * (param_mr_str.tau_b + param_user_str.delta_bolus);
			bolus_time_passed = bolus_arrived * delta_bolus;

			delta_M_tissue = delta_M_tissue + calculate_M0_tissue_Hrabe_no_dispersion(cbf, tau, t, t1_tissue, t1_blood, arrival_time, bolus_time_passed);

			bolus_arrived = bolus_arrived + 1;
		end
	end


end

