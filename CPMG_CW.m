% GNU Octave / MATLab script to simulate relaxation dispersion experiment
% SQ in-phase N15 CPMG using matched CW 1H decoupling (Flemmming's experiment)

% HvI 2016 from CPMG_CW

if feedback >= 1
    disp("")
    disp(" *** CW-matched CPMG parameters ***")
end

% pre-define some propagators
pNx90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi/(2*gB1_N));   % 90 x on N, used to excite magnetization
pNx90m    = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi/(2*gB1_N));   % 90 -x on N

pNx180    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi/(gB1_N));     % 180 x on N, used in P-element
pNx180m   = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi/(gB1_N));     % 180 -x on N

pHx180    = expm(LVM(offset, pi*JNH, gB1_H*[1 0  0 0], relaxRates)*pi/(gB1_H));     % 180 x on H, used in P-element

%============================================================================
%
% now do simulation
%
%============================================================================

for kk=0:npoints                    % 0=reference experiment, >=1 is dispersion experiment


    if kk==0                        % do reference experiment
        cv = StartOp;               % start from Nz
        cv = pNx90p*cv;             % only 90x pulse
        if cycle_phi2 == 1
            cv_scan1 = cv;                  % store scan1
            cv_scan2 = pNx90m*StartOp;      % only 90p
            cv       = 0.5*(cv_scan1 - cv_scan2); % apply phase cycle
        end
    end

    if (kk>0)               % do 90 w/ CW and first half CT on -Ny

        tau  = 1/(4*nuCPMG(kk)) - pw_cpmg*1e-6;     % since we now simulate w/ finite pulse length need to compensate for pulse length!
        tau1 = tau - 0.636*pw_90*1e-6;              % compensate first tau_cp for evolution during excite pulse to align magnetization
        k = numPulses(kk)/2;                        % numPulses is total number of pulses, k is number of pulses per CPMG half
 
        step    = round(round(0.5/(pw_dec*4)*1e4)*1e2/nuCPMG(kk));
        gB1_cw  = 2*step*nuCPMG(kk)*2*pi;               % set wB1 for 1H CW decoupling; should be roughly 15 kHz / 16-17 us
        if feedback == 2
            printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_cw/(2*pi), pi/(2*gB1_cw)*1e6)
        elseif feedback == 1 && kk==1
            printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_cw/(2*pi), pi/(2*gB1_cw)*1e6)
        end

        pNy180_cw = expm(LVM(offset, pi*JNH, gB1_cw*[1 0 0 0] + gB1_cp*[0 0  0 1], relaxRates)*pi/(gB1_cp));   % 180 Ny pulse + 1H CW decoupling!
        LVdecexp  = expm(LVM(offset, pi*JNH, gB1_cw*[1 0 0 0], relaxRates)*tau);                    % store the matrix exponential for increased speed
        LVdecexp1 = expm(LVM(offset, pi*JNH, gB1_cw*[1 0 0 0], relaxRates)*tau1);                   % store the matrix exponential for increased speed

        % 1. excite Nz magnetization and start CW at evolution point of 90
        cv = StartOp;               % start from Nz
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi*0.363/(2*gB1_N))*cv;                     % first apply 15N pulse of length 0.363*pw_90
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0] + gB1_cw*[1 0 0 0], relaxRates)*pi*0.636/(2*gB1_N))*cv;  % remainder of pulse with 1H dec
        if cycle_phi2 == 1
            cv_scan1 = cv;                                                                                      % store scan1
            cv_scan2 = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi*0.363/(2*gB1_N))*StartOp;      % first apply 15N pulse of length 0.363*pw_90
            cv_scan2 = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0] + gB1_cw*[1 0 0 0], relaxRates)*pi*0.636/(2*gB1_N))*cv_scan2;   % remainder of pulse with 1H dec
            cv       = 0.5*(cv_scan1 - cv_scan2);                                                               % apply phase cycle
        end
        cv_start = cv;                                                                                          % store starting point for different nuCPMG values

        % 2.now evolve magnetization during tau_cp -180 - tau_cp spin echos
        cv = LVdecexp1*cv_start;                        % first tau_cp is compensated for evolution in pulse
        cv = pNy180_cw*cv;                              % apply 180 along y
        cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
        for l=1:k-1                                     % for remaining (tau-180-tau) element do:
            cv = LVdecexp*cv;                           % free precession with exchange during tau
            cv = pNy180_cw*cv;                          % apply 180 along y
            cv = LVdecexp*cv;                           % free precession with exchange during tau
        end                                             % iteration over k

    end                     % end first half CT

    % now apply P-element, also in reference = delay- 180N,Hx-delay

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
        %cv = cv_before; % to check influence of P-element on slow pulsing artifacts -- YES this is the critical factor., in fact it is the 180 pulse!
        for l=1:k-1                                     % for remaining (tau-180-tau) element do:
            cv = LVdecexp*cv;                           % evolve tau with exchange
            cv = pNy180_cw*cv;                          % apply 180 along y
            cv = LVdecexp*cv;                           % evolve tau with exchange
        end                                             % iteration over k
        cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
        cv = pNy180_cw*cv;                              % apply 180 along y
        cv = LVdecexp1*cv;                              % first tau_cp is compensated for evolution in pulse

        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0]+gB1_cw*[1 0 0 0], relaxRates)*pi*0.636/(2*gB1_N))*cv;    % remainder of pulse with 1H dec
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi*0.363/(2*gB1_N))*cv;                     % 90 -x on N, 0.363*p21
    end                     % end second half CT

    if kk==0                    % reference experiment
        cv = pNx90m*cv;         % return to z with 90mx (to match catia)
    end

    LV_eq = expm(LVM(offset, pi*JNH, 0*[1 0 0 0], relaxRates)*tau_eq);      % equilibrate
    cz = LV_eq*cv;          

    % calculate observables
    mag(kk+1) = real(cz(7));        % Nz (there is no ground/excited state in this setup)
    if (kk>0)
        R2eff(kk) = -(1/time_T2)*log(mag(kk+1)/mag(1));
    end

end	%kk

R2_CW = R2eff;

