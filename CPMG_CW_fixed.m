% GNU Octave / MATLab script to simulate strong H-H coupling artefacts in 
% TROSY/anti-TROSY CPMG experiments
%
% Hvi 2016 from CPMG_CW.m   
%           now using a single, fixed CW power for all nuCPMG values

if feedback >= 1
    disp("")
    disp(" *** CW-fixed CPMG parameters ***")
end

% pre-define some propagators
pNx90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi/(2*gB1_N));   % 90 x on N
pNx90m    = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi/(2*gB1_N));   % 90 -x on N
pNy90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  0 1], relaxRates)*pi/(2*gB1_N));   % 90 y on N

pNx180    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi/(gB1_N));     % 180 x on N
pNx180m   = expm(LVM(offset, pi*JNH, gB1_N*[0 0  -1 0], relaxRates)*pi/(gB1_N));    % 180 -x on N

pHx180    = expm(LVM(offset, pi*JNH, gB1_H*[1 0  0 0], relaxRates)*pi/(gB1_H));     % 180 x on H

pNy180_cw = expm(LVM(offset, pi*JNH, gB1_cp*[0 0  0 1]+gB1_H*[1 0 0 0], relaxRates)*pi/(gB1_cp));   % 180 Ny pulse + 1H CW decoupling!

%============================================================================
%
% now do simulation
%
%============================================================================

cv = StartOp;               % define density matrix for ground and excited state
cv = pNx90p*cv;             % aply 90x on Nz --- need to check whether phase cycle is necessary
cv_start = cv;              % store starting point for different nuCPMG values

if cycle_phi2 == 1
	cv_scan1 = cv;                  % store scan1
	cv_scan2 = pNx90m*StartOp;      % aply 90-x on Nz in scan2
	cv       = 0.5*(cv_scan1 - cv_scan2); % apply phase cycle
	cv_start = cv;                  % overwrite cv_start
end

for kk=0:npoints            % 0=reference experiment, >=1 is dispersion experiment

	if (kk>0)               % do first half CT on -Ny

		if offResonance == 1
			tau  = 1/(4*nuCPMG(kk)) - pw_cpmg*1e-6;     % since we now simulate w/ finite pulse length need to compensate for pulse length!
			tau1 = tau - 0.636*pw_90*1e-6;              % compensate first tau_cp for evolution during excite pulse to align magnetization
		else
			tau = 1/(4*nuCPMG(kk));                     % set tauCPMG 2*k*2*tau=CT <=> tau = CT/(4k)
			tau1= tau;
		end
		k = numPulses(kk)/2;                            % numPulses is total number of pulses

		if feedback == 2
			printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_H/(2*pi), pi/(2*gB1_H)*1e6)
		elseif feedback == 1 && kk==1
			printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_H/(2*pi), pi/(2*gB1_H)*1e6)
		end
		LVdecexp  = expm(LVM(offset, pi*JNH, gB1_H*[1 0 0 0], relaxRates)*tau);                    % store the matrix exponential for increased speed
		LVdecexp1 = expm(LVM(offset, pi*JNH, gB1_H*[1 0 0 0], relaxRates)*tau1);                   % store the matrix exponential for increased speed

		% now evolve magnetization
		cv = LVdecexp1*cv_start;                        % first tau_cp is compensated for evolution in pulse
		cv = pNy180_cw*cv;                              % apply 180 along y
		cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
		for l=1:k-1                                     % for remaining (tau-180-tau) element do:
			cv = LVdecexp*cv;                           % free precession with exchange during tau
			cv = pNy180_cw*cv;                          % apply 180 along y
			cv = LVdecexp*cv;                           % free precession with exchange during tau
		end                                             % iteration over k

	end                     % end first half CT

	% now apply P-element, also in reference = - 180N,Hx-

	cv_before = cv;
	cv = pNx180*cv;         % apply 180 N and H
	cv = pHx180*cv;

	if cycle_phi3 == 1                                  % Jiang et al. use phase cycle here!, Flemming states x,-x
		cv_scan1 = cv;                                  % store scan1
		cv       = pNx180m*cv_before;                   % do scan 2, apply 180 N and H
		cv_scan2 = pHx180*cv;
		cv       = 0.5*(cv_scan1 + cv_scan2);           % do cycling
	end

	% end P-element

	if (kk>0)               % now apply second half of CT-CPMG again on Ny

		for l=1:k-1                                     % for remaining (tau-180-tau) element do:
			cv = LVdecexp*cv;                           % evolve tau with exchange
			cv = pNy180_cw*cv;                          % apply 180 along y
			cv = LVdecexp*cv;                           % evolve tau with exchange
		end                                             % iteration over k
		cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
		cv = pNy180_cw*cv;                              % apply 180 along y
		cv = LVdecexp1*cv;                              % first tau_cp is compensated for evolution in pulse

	end                     % end second half CT

	cz = pNx90m*cv;         % return to z with 90mx (to match catia)
	LV_eq = expm(LVM(offset, pi*JNH, 0*[1 0 0 0], relaxRates)*tau_eq);
	cz = LV_eq*cz;          % equilibrate
	
	% calculate observables; t	% calculate observables

	mag(kk+1) = real(cz(7));        % Nz major


  	REF = mag(1);
   	if (kk>0)
		R2eff(kk) = -(1/time_T2)*log(mag(kk+1)/REF);
   	end

end	%kk

R2_CWf = R2eff;
