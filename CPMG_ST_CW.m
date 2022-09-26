% GNU Octave / MATLab script to simulate relaxation dispersion experiment
% SQ in-phase N15 CPMG using fixed CW 1H decoupling and Zuiderweg phase cycle (Jiang's experiment)

% HvI 2016

if feedback >= 1
    disp("")
    disp(" *** ST-CW CPMG parameters ***")
end

% pre-define some propagators
pNx90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*pi/(2*gB1_N));   % 90 +x on N
pNx90m    = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi/(2*gB1_N));   % 90 -x on N
pNy90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  0 1], relaxRates)*pi/(2*gB1_N));   % 90 +y on N

pNy180p_cw = expm(LVM(offset, pi*JNH, gB1_cp*[0 0  0 1]+gB1_H*[1 0 0 0], relaxRates)*pi/(gB1_cp));   % 180 N+y pulse + 1H CW decoupling!
pNx180p_cw = expm(LVM(offset, pi*JNH, gB1_cp*[0 0  1 0]+gB1_H*[1 0 0 0], relaxRates)*pi/(gB1_cp));   % 180 N+x pulse + 1H CW decoupling!
pNx180m_cw = expm(LVM(offset, pi*JNH, gB1_cp*[0 0 -1 0]+gB1_H*[1 0 0 0], relaxRates)*pi/(gB1_cp));   % 180 N-x pulse + 1H CW decoupling!

Zuiderweg  = zeros(16,16,4);
Zuiderweg(:,:,1)  = pNy180p_cw;
Zuiderweg(:,:,2)  = pNy180p_cw;
Zuiderweg(:,:,3)  = pNx180p_cw;
Zuiderweg(:,:,4)  = pNx180m_cw;                                      % define 3D array corresponding to Zuiderweg phase cycle
Zuidertext = cellstr(["+y";"+y";"+x";"-x"]);

% test effect of XY phasecycle on slow pulsing artifact
%Zuidertext = cellstr(["+y";"+y";"+y";"+y"]);
%Zuiderweg(:,:,3)  = pNy180p_cw;
%Zuiderweg(:,:,4)  = pNy180p_cw;    
%NO EFFECT, slow pulsing artifact is still larger than in Flemming's version

pHx180    = expm(LVM(offset, pi*JNH, gB1_H*[1 0  0 0], relaxRates)*pi/(gB1_H));     % 180 x on H, used in spin-echo on reference

gB1_cw = gB1_H;

%============================================================================
%
% now do simulation
%
%============================================================================

for kk=0:npoints            % 0=reference experiment, >=1 is dispersion experiment

    if kk==0                        % do reference experiment
        cv = StartOp;               % start from Nz
        cv = pNx90p*cv;             % only 90x pulse
        if cycle_phi2 == 1
            cv_scan1 = cv;                  % store scan1
            cv_scan2 = pNx90m*StartOp;      % only 90p
            cv       = 0.5*(cv_scan1 - cv_scan2); % apply phase cycle
        end

        LV_4us = expm(LVM(offset, pi*JNH, 0*[1 0 0 0], relaxRates)*4e-6);
        cv = LV_4us*cv;
        cv = pHx180*cv;
        cv = pNy90p*cv;             % 180 at high power
        cv = pNy90p*cv;
        cv = LV_4us*cv;

        cv = pNx90m*cv;             % return to z with 90mx (to match catia)
    end

    if (kk>0)               % do single-train CPMG

        tau  = 1/(4*nuCPMG(kk)) - pw_cpmg*1e-6;     % since we now simulate w/ finite pulse length need to compensate for pulse length!
        tau1 = tau - 0.636*pw_90*1e-6;              % compensate first tau_cp for evolution during excite pulse to align magnetization
        if feedback == 2
            printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_H/(2*pi), pi/(2*gB1_H)*1e6)
        elseif feedback == 1 && kk==1
            printf("\tvCPMG %6.1f, tau_cp %6.3f us, dec. @ %6.1f kHz (%5.2f us 90)\n", nuCPMG(kk), tau*1e6, gB1_H/(2*pi), pi/(2*gB1_H)*1e6)
        end
        LVdecexp  = expm(LVM(offset, pi*JNH, gB1_H*[1 0 0 0], relaxRates)*tau);                    % store the matrix exponential for increased speed
        LVdecexp1 = expm(LVM(offset, pi*JNH, gB1_H*[1 0 0 0], relaxRates)*tau1);                   % store the matrix exponential for increased speed

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

        % 2. now evolve magnetization during tau_cp -180 - tau_cp spin echos, coded as in pulse sequence
        if numPulses(kk) == 2
            cv = LVdecexp1*cv_start;                        % first tau_cp is compensated for evolution in pulse
            cv = pNy180p_cw*cv;                             % apply 180 along y, this is actually phase_cycled y,x
            cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
            cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
            cv = pNy180p_cw*cv;                             % apply 180 along y, this is actually phase_cycled y,-x
            cv = LVdecexp1*cv;                              % last tau_cp is compensated for evolution in pulse
            if cycle_phi7 == 1
                cv_scan1 = cv;
                cv = LVdecexp1*cv_start;                    % first tau_cp is compensated for evolution in pulse
                cv = pNx180p_cw*cv;                         % apply 180 along y, this is actually phase_cycled y,x
                cv = LVdecexp*cv;                           % free precession with exchange during tau_cp
                cv = LVdecexp*cv;                           % free precession with exchange during tau_cp
                cv = pNx180m_cw*cv;                         % apply 180 along y, this is actually phase_cycled y,-x
                cv = LVdecexp1*cv;                          % last tau_cp is compensated for evolution in pulse
                cv = 0.5*(cv + cv_scan1);                   % apply phase-cycle
            end
        elseif numPulses(kk) > 2
            ph19 = 1;                                       % reset phase counter
            cv = LVdecexp1*cv_start;                        % first tau_cp is compensated for evolution in pulse
            cv = Zuiderweg(:,:,mymod(ph19,4))*cv;           % apply 180 along y, this is actually phase_cycled y,x
            cv = LVdecexp*cv;                               % free precession with exchange during tau_cp
            ph19 += 1;                                      % increment phase

            k  = 0.5*(numPulses(kk)-2);                     % first & last pulse is explicit, two pulses per loop
            for l=1:k                                       % for remaining (tau-180-tau) element do:
                cv = LVdecexp*cv;                           % free precession with exchange during tau
                cv = Zuiderweg(:,:,mymod(ph19,4))*cv;       % apply 180 along y
                cv = LVdecexp*cv;                           % free precession with exchange during tau
                ph19 += 1;                                  % increment phase
                cv = LVdecexp*cv;                           % free precession with exchange during tau
                cv = Zuiderweg(:,:,mymod(ph19,4))*cv;       % apply 180 along y
                cv = LVdecexp*cv;                           % free precession with exchange during tau
                ph19 += 1;                                  % increment phase
            end                                             % iteration over k

            cv = LVdecexp*cv;                               % free precession with exchange during tau
            cv = Zuiderweg(:,:,mymod(ph19,4))*cv;           % apply 180 along y
            cv = LVdecexp1*cv;                              % free precession with exchange during tau compensated
        end  % end of numpulses >2
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0]+gB1_cw*[1 0 0 0], relaxRates)*pi*0.636/(2*gB1_N))*cv;    % remainder of pulse with 1H dec
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*pi*0.363/(2*gB1_N))*cv;                     % 90 -x on N, 0.363*p21
    end % end dispersion plane, do reference

    % equilibrate
    LV_eq = expm(LVM(offset, pi*JNH, 0*[1 0 0 0], relaxRates)*tau_eq);
    cz = LV_eq*cv;          % equilibrate

    % calculate observables
    mag(kk+1) = real(cz(7));        % Nz (there is no ground/excited state in this setup)
    if (kk>0)
        R2eff(kk) = -(1/time_T2)*log(mag(kk+1)/mag(1));
    end

end	%kk

R2_ST  = R2eff;

