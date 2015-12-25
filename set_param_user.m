% This code is the implementation the following papers
% MA Chappell (2012) doi: 10.1002/mrm.24372 (MACQ) QUASAR paper
% MA Chappell (2012) doi: 10.1002/mrm.24260 (MACD) Dispersion paper
% MA Chappell (2011) doi: 10.1002/mrm.22641 (MACP) Partial volume paper
% ET Petersen (2006) doi: 10.1002/mrm.20784 (ETP)

function [] = set_param_user(varargin)

	% Default variables
	% User input variables
	param_user_str                       = struct;
	param_user_str.f                     = 100 / 6000; % CBF need to convert from ml/100ml/min to ml/g/s
	param_user_str.f_gm                  = 60 / 6000; % GM CBF needed to convert from ml/100ml/min to ml/g/s
	param_user_str.f_wm                  = 20 / 6000; % WM CBF needed to convert from ml/100ml/min to ml/g/s
	param_user_str.arterial_blood_volume = 1.5 / 100; % arterial blood volume in percentage
	param_user_str.tau_t                 = 0; % Bolus arrival time to tissue (microvasculature)
	param_user_str.tau_t_gm              = 0.5; % Bolus arrival time to GM tissue
	param_user_str.tau_t_wm              = 1.0; % Bolus arrival time to WM tissue
	param_user_str.tau_m                 = 0.5; % Bolus arrival time to vasculature
	param_user_str.t1_t                  = 1.3; % T1 relaxation of tissue
	param_user_str.t1_t_gm               = 1.3; % T1 relaxation time of GM tissue
	param_user_str.t1_t_wm               = 1.1; % T1 relaxation time of WM tissue
	param_user_str.t1_a                  = 1.6; % T1 relaxation of arterial blood, in postprocessing section of (ETP)
	param_user_str.t1_t_gm_correct       = 0.68; % Look-locker corrected T1 relaxation of GM tissue
	param_user_str.t1_t_wm_correct       = 0.62; % Look-locker corrected T1 relaxation of WM tissue
	param_user_str.t1_a_correct          = 0.76; % Look-locker corrected T1 relaxation of arterial blood
	param_user_str.delta_ti              = 0.6; % interval between excitation pulses (this include the previous bolus duration)
	param_user_str.delta_ti_gap_factor   = 1; % number of delta_TIs between each successive bolus (default is 1)
	param_user_str.delta_bolus           = param_user_str.delta_ti_gap_factor * param_user_str.delta_ti; % temporal duration between successive boluses (this include the previous bolus duration)
	param_user_str.slice_shifting_factor = 2; % Factors of slice shifting
	param_user_str.actual_sampling_rate  = param_user_str.delta_ti / param_user_str.slice_shifting_factor;
	param_user_str.t                     = 0.04 : param_user_str.actual_sampling_rate : 7; % sampling time points, second variable must be equal to actual sampling rate (Default: 11 TIs)
	param_user_str.m_0a                  = 1; % equilibrium magnetization of arterial blood
	param_user_str.inversion_efficiency  = 0.91; % inversion efficiency alpha, in postprocessing section of (ETP)
	param_user_str.snr                   = 40; % signal to noise ratio (in percentage)
	param_user_str.sd                    = 1000; % standard deviation of noise
	param_user_str.crush_efficiency      = 0; % percentage of arterial blood signal remaining after applying crusher gradients
	
	param_user_str.dispersion_type       = 1; % no dispersion
	
	param_user_str.mask                  = 'ref_images/mask_mean'; % file name of Mask file
	param_user_str.pvgm                  = 'ref_images/pvgm_one'; % 'mask_bin'; % file name of PV GM map
	param_user_str.pvwm                  = 'ref_images/pvwm_zero'; % 'mask_bin'; % file name of PV WM map
	param_user_str.abv_mask              = 'ref_images/mask_abv'; % 'mask_abv'; file name of ABV mask (currently all zeros)

	if(length(varargin) == 1)
		current_handles                      = varargin{1};
		param_user_str.f                     = str2double(get(current_handles.et_cbf, 'String')) / 6000;
		param_user_str.arterial_blood_volume = str2double(get(current_handles.et_abv, 'String')) / 100;
		param_user_str.tau_t                 = str2double(get(current_handles.et_tau_t, 'String'));
		param_user_str.tau_m                 = str2double(get(current_handles.et_tau_m, 'String'));
		param_user_str.t1_t                  = str2double(get(current_handles.et_t1_t, 'String'));
		param_user_str.t1_a                  = str2double(get(current_handles.et_t1_a, 'String'));
		param_user_str.t                     = eval(get(current_handles.et_ti, 'String'));

	end


	% save the user parameters in a file
	filename = 'param_user.mat';
	save(filename, 'param_user_str');

end

