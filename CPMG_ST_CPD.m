% GNU Octave / MATLab script to simulate relaxation dispersion experiment
% SQ in-phase N15 CPMG using matched CPD 1H decoupling and Zuiderweg phase cycle (our experiment)
%
% HvI Nov 2016

if feedback >= 1
    disp("")
    printf(" *** matched ST-CPD (using %s CPD) CPMG parameters ***\n", cpdString{cc+1})
end

% calculate pulse lengths 15N
tp_N    = pi/(2*gB1_N);                 % = pw_90*1e-6   ~37us
tp_cpmg = pi/(gB1_cp);                  % = pw_cpmg*1e-6 ~45us

% derive pulse lengths CPD element 1/2/3
tp_1 = cpd_angle(1)/90*pi/(2*gB1_H);  % ~17us
tp_2 = cpd_angle(2)/90*pi/(2*gB1_H);  % ~45us ==>. CPD element = 80 us., each element is less than a 15N 180
tp_3 = cpd_angle(3)/90*pi/(2*gB1_H);  % ~17us
tp_CPD  = tp_1 + tp_2 + tp_3;

% pre-define some propagators
% excitation pulse 15N
pNx90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  1 0], relaxRates)*tp_N);   % 90 x on N
pNx90m    = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0], relaxRates)*tp_N);   % 90 -x on N
pNy90p    = expm(LVM(offset, pi*JNH, gB1_N*[0 0  0 1], relaxRates)*tp_N);   % 90 y on N

% define LV-container w/ input for LV_CPD function, these are all constants within one simulation
LVpars{1} = offset;
LVpars{2} = pi*JNH;
LVpars{3} = relaxRates;
LVpars{4} = cpd_angle;
global remTime

% find the matched 1H decoupling power such that each tau-90 has an integer number of CPD blocks
% otherwise there will be huge spikes in the simulations (not really seen experimentally)
tauValues         = 1./(4*nuCPMG);
numberOfBlocks    = floor(tauValues./tp_CPD)+1;
cpdLength         = tauValues./numberOfBlocks;
matchedDec        = 1e6*cpdLength/(sum(cpd_angle)/90);
newBlockFractions = (1./(4*nuCPMG) )./cpdLength - floor((1./(4*nuCPMG) )./cpdLength);

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
            cv_scan1 = cv;                        % store scan1
            cv_scan2 = pNx90m*StartOp;            % only 90-x
            cv       = 0.5*(cv_scan1 - cv_scan2); % apply phase cycle
        end

        LV_4us = expm(LVM(offset, pi*JNH, 0*[1 0 0 0], relaxRates)*4e-6);
        cv = LV_4us*cv;
        cv = pHx180*cv;             % spin echo on water 
        cv = pNy90p*cv;             % 180 at high power 
        cv = pNy90p*cv;
        cv = LV_4us*cv;

        cv = pNx90m*cv;             % return to z with 90mx (to match catia)
    end

    if (kk>0)               % do single-train CPMG

        tau    = 1/(4*nuCPMG(kk)) - pw_cpmg*1e-6;     % since we now simulate w/ finite pulse length need to compensate for pulse length!
        tau1   = tau - 0.636*pw_90*1e-6;              % compensate first tau_cp for evolution during excite pulse to align magnetization
        gB1_cw = 2*pi*1/(4e-6*matchedDec(kk) ) ;      % set wB1 for to matched 1H CPD decoupling power

        % derive pulse lengths CPD element 1/2/3  --- needs to done for each nuCPMG value!
        global tp_1 tp_2 tp_3 tp_CPD
        tp_1 = cpd_angle(1)/90*pi/(2*gB1_cw);  % ~17us
        tp_2 = cpd_angle(2)/90*pi/(2*gB1_cw);  % ~45us ==>. CPD element = 80 us., each element is less than a 15N 180
        tp_3 = cpd_angle(3)/90*pi/(2*gB1_cw);  % ~17us
        tp_CPD  = tp_1 + tp_2 + tp_3;

        % pre-define 9 [V,D]-pairs for partial CPD elements for speed and make them globally accessible
        global V1_0 V2_0 V3_0 V1_Nyp V2_Nyp V3_Nyp V1_Nxp V2_Nxp V3_Nxp V1_Nxm V2_Nxm V3_Nxm
        global D1_0 D2_0 D3_0 D1_Nyp D2_Nyp D3_Nyp D1_Nxp D2_Nxp D3_Nxp D1_Nxm D2_Nxm D3_Nxm

        [V1_0  ,D1_0  ] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(1,:), relaxRates));                            %  90 x on H, no pulse on 15N
        [V2_0  ,D2_0  ] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(2,:), relaxRates));                            % 240 y on H, no pulse on 15N
        [V3_0  ,D3_0  ] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(3,:), relaxRates));                            %  90 x on H, no pulse on 15N

        [V1_Nyp,D1_Nyp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(1,:) + gB1_cp*[0 0  0 1], relaxRates));        %  90 x on H, 15N y
        [V2_Nyp,D2_Nyp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(2,:) + gB1_cp*[0 0  0 1], relaxRates));        % 240 y on H, 15N y
        [V3_Nyp,D3_Nyp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(3,:) + gB1_cp*[0 0  0 1], relaxRates));        %  90 x on H, 15N y

        [V1_Nxp,D1_Nxp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(1,:) + gB1_cp*[0 0  1 0], relaxRates));        %  90 x on H, 15N x
        [V2_Nxp,D2_Nxp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(2,:) + gB1_cp*[0 0  1 0], relaxRates));        % 240 y on H, 15N x
        [V3_Nxp,D3_Nxp] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(3,:) + gB1_cp*[0 0  1 0], relaxRates));        %  90 x on H, 15N x

        [V1_Nxm,D1_Nxm] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(1,:) + gB1_cp*[0 0 -1 0], relaxRates));        %  90 x on H, 15N -x
        [V2_Nxm,D2_Nxm] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(2,:) + gB1_cp*[0 0 -1 0], relaxRates));        % 240 y on H, 15N -x
        [V3_Nxm,D3_Nxm] = eig(LVM(offset, pi*JNH, gB1_cw*cpdPhase(3,:) + gB1_cp*[0 0 -1 0], relaxRates));        %  90 x on H, 15N -x
        
        % pre-define 9 full CPD elements for speed and make them globally accessible
        global CPD1 CPD2 CPD3 CPD1_Nyp CPD2_Nyp CPD3_Nyp CPD1_Nxp CPD2_Nxp CPD3_Nxp CPD1_Nxm CPD2_Nxm CPD3_Nxm
        
        CPD1      = real(V1_0*diag(exp(diag(D1_0)*tp_1))/V1_0);                            %  90 x on H, no pulse on 15N
        CPD2      = real(V2_0*diag(exp(diag(D2_0)*tp_2))/V2_0);                            % 240 y on H, no pulse on 15N
        CPD3      = real(V3_0*diag(exp(diag(D3_0)*tp_3))/V3_0);                           %  90 x on H, no pulse on 15N

        CPD1_Nyp  = real(V1_Nyp*diag(exp(diag(D1_Nyp)*tp_1))/V1_Nyp);                      %  90 x on H, 15N y
        CPD2_Nyp  = real(V2_Nyp*diag(exp(diag(D2_Nyp)*tp_2))/V2_Nyp);;                     % 240 y on H, 15N y
        CPD3_Nyp  = real(V3_Nyp*diag(exp(diag(D3_Nyp)*tp_3))/V3_Nyp);;                     %  90 x on H, 15N y

        CPD1_Nxp  = real(V1_Nxp*diag(exp(diag(D1_Nxp)*tp_1))/V1_Nxp);                      %  90 x on H, 15N x
        CPD2_Nxp  = real(V2_Nxp*diag(exp(diag(D2_Nxp)*tp_2))/V2_Nxp);                      % 240 y on H, 15N x
        CPD3_Nxp  = real(V3_Nxp*diag(exp(diag(D3_Nxp)*tp_3))/V3_Nxp);                      %  90 x on H, 15N x

        CPD1_Nxm  = real(V1_Nxm*diag(exp(diag(D1_Nxm)*tp_1))/V1_Nxm);                      %  90 x on H, 15N -x
        CPD2_Nxm  = real(V2_Nxm*diag(exp(diag(D2_Nxm)*tp_2))/V2_Nxm);                      % 240 y on H, 15N -x
        CPD3_Nxm  = real(V3_Nxm*diag(exp(diag(D3_Nxm)*tp_3))/V3_Nxm);                      %  90 x on H, 15N -x

        % pre-define 4 full CPD block propagators for speed and make them globally accessible
        global CPDt CPDt_Nyp CPDt_Nxp CPDt_Nxm
        CPDt     = CPD3    *CPD2    *CPD1;
        CPDt_Nyp = CPD3_Nyp*CPD2_Nyp*CPD1_Nyp;
        CPDt_Nxp = CPD3_Nxp*CPD2_Nxp*CPD1_Nxp;
        CPDt_Nxm = CPD3_Nxm*CPD2_Nxm*CPD1_Nxm;

        if feedback == 2
            printf("\tvCPMG %6.1f, tau-90 %6.3f us, tp_CPD @ %.3f us, (%5.3f blocks)\n", nuCPMG(kk), 1e6/(4*nuCPMG(kk)), 1e6*tp_CPD, (1/(4*nuCPMG(kk)))/tp_CPD)
        elseif feedback == 1 && kk==1
            printf("\tvCPMG %6.1f, tau-90 %6.3f us, tp_CPD @ %.3f us, (%5.3f blocks)\n", nuCPMG(kk), 1e6/(4*nuCPMG(kk)), 1e6*tp_CPD, (1/(4*nuCPMG(kk)))/tp_CPD)
        end

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
            cv = LVM_CPD(cv_start,         tau1,       0,                0, LVpars);             % first tau_cp is compensated for evolution in pulse
            cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*[0 0 0 1], LVpars);             % apply 180 along y, this is actually phase_cycled y,x
            cv = LVM_CPD(cv      ,         2*tau, remTime,               0, LVpars);             % free precession with exchange during 2*tau_cp
            cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*[0 0 0 1], LVpars);             % apply 180 along y, this is actually phase_cycled y,x
            cv = LVM_CPD(cv      ,         tau1, remTime,                0, LVpars);             % last tau_cp is compensated for evolution in pulse
            if cycle_phi7 == 1
                cv_scan1 = cv;
                cv = LVM_CPD(cv_start,         tau1,       0,                 0, LVpars);        % first tau_cp is compensated for evolution in pulse
                cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*[0 0  1 0], LVpars);        % apply 180 along y, this is actually phase_cycled y,x
                cv = LVM_CPD(cv      ,        2*tau, remTime,                 0, LVpars);        % free precession with exchange during 2*tau_cp
                cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*[0 0 -1 0], LVpars);        % apply 180 along y, this is actually phase_cycled y,x
                cv = LVM_CPD(cv      ,         tau1, remTime,                 0, LVpars);        % last tau_cp is compensated for evolution in pulse
                cv = 0.5*(cv + cv_scan1);                                                        % apply phase-cycle
            end
        elseif numPulses(kk) > 2
            ph19 = 1;                                                                            % reset phase counter
            cv = LVM_CPD(cv_start,         tau1,       0,                      0, LVpars);       % first tau_cp is compensated for evolution in pulse
            cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*Zuidersel(ph19), LVpars);       % apply 180 along y
            cv = LVM_CPD(cv      ,        tau, remTime,                      0, LVpars);         % free precession with exchange during tau_cp
            ph19 += 1;                                                                           % increment phase
            k  = (numPulses(kk)-2);                                                              % first & last pulse is explicit, one pulses per loop
            for l=1:k                                                                            % for remaining (tau-180-tau) element do:
                cv  = LVM_CPD(cv      ,        tau, remTime,                      0, LVpars);    % free precession with exchange during tau_cp
                cv  = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*Zuidersel(ph19), LVpars);  % apply 180 along Zuiderweg cycle
                cv  = LVM_CPD(cv      ,        tau, remTime,                      0, LVpars);    % free precession with exchange during tau_cp
                ph19 = ph19+1;                                                                   % increment phase
            end                                         % iteration over k
            cv = LVM_CPD(cv      ,        tau, remTime,                      0, LVpars);         % free precession with exchange during tau_cp
            cv = LVM_CPD(cv      , pw_cpmg*2e-6, remTime, gB1_cp*Zuidersel(ph19), LVpars);       % apply 180 along Zuiderweg cyc
            cv = LVM_CPD(cv      ,         tau1, remTime,                      0, LVpars);       % last tau_cp is compensated for evolution in pulse
        end % numPulses
        cv = expm(LVM(offset, pi*JNH, gB1_N*[0 0 -1 0]+gB1_cw*[1 0 0 0], relaxRates)*pi*0.636/(2*gB1_N))*cv;   % remainder of pulse with 1H dec
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

R2_CPDm = R2eff;
