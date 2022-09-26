% GNU Octave / MATLab script to simulate 15N CPMG
%
% definePars.m --- Allard style, no exchange!
% 
% define basic parameters

offResonance        = 1;    % [0|1] include offResonance yes/no
includeRelaxation   = 1;    % [0|1] include relaxation yes/no
plot_hold           = 0;    % [0|1] delete | keep old plots -- with automatic color switching
plot_jiang          = 0;    % [0|1] make plot as in Fig 3 and S4 of Jiang et al.
feedback            = 1;    % [0|1|2] no | yes | all give some feedback on parameters
cycle_phi3          = 1;    % [0|1] do phase cycling on central 180 15N no (0) yes (1)
cycle_phi2          = 1;    % [0|1] do phase cycling on excite 90 15N no (0) yes (1)
cycle_phi7          = 1;    % [0|1] do phase cycling on ST-CW numPulses =2 CPMG pulse no (0) yes (1)

% spin system
JNH     = -92;                             % 1JNH [Hz]
wNppm   = 0/81.1;                          % 15N offset in ppm wrt to 15N carrier
wHppm   = 3;                               % 1H offset in ppm wrt to 1H CW carrier

% experimental set-up
pw_cpmg = 45;                                                   % 15N 90° pulse length CPMG
pw_90   = 30;                                                   % 15N 90 high power for excite Nz to Ny and for bringing back to Nz
pw_dec  = 17;                                                   % 1H 90° decoupling pulse, 
cfH     = 1;                                                    % calibration factor 1H pulses: 1 = perfectly calibrated, 1.1 = +10% 0.9=-10%
cfN     = 180/180;                                              % calibration factor 15N pulses: 1 = perfectly calibrated, 1.1 = +10% 0.9=-10%
B0      = 19.964;                                               % static magnetic field corresponding to 800 MHz [T]
                                                                % 500: 11.7434 600: 14.0921 750: 800: 18.789 850: 19.964 950: 1.2GHz:
time_T2 = 0.04;                                                 % length of CT CPMG relaxation period [s]
maxCPMG = 1000;                                                 % maximum CPMG field [Hz]
tau_eq  = 0.003e-3;	                                            % equilibration time at end of sequence 

% define dynamics
tc = 9e-9;            % rotational correlation time [s]
                      % Jiang et al. uses 9ns
S2 = 0.85;            % generalized order parameter squared
                      % not specified in Jiang et al. but 0.85 coves 15N R2 of 13.7 as in their case
te = 0.1e-9;          % correlation time internal motions [s] (fixed at 100 ps)

% output files
outfile = 'test';

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%            DO NOT CHANGE BELOW THIS LINE
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% fixed parameters
gH      = 2.6752e8;         % 1H gyromagnetic ratio [ T-1 s-1]
gN      = -2.7120e7;        % 15N gyromagnetic ratio [ T-1 s-1]
hbar    = 1.0546e-34;       % Planck's constant over 2*pi [J s ]
CSAN    = -172e-6;          % 15N chemical shift anisotropy [ppm]
CSAH    = 10e-6;            % 1H  chemical shift anisotropy [ppm]  (10 ppm seems to be bit large just as effective as 15N ccr?)
rNH     = 1.02e-10;         % distance N-H [A]
rHH     = 1.85e-10;         % 1HN - 1Halpha distance [A] (effective value see Jiang paper)
                            % 2.25 A corresponds to 13.7 / 17.1 inphase/antiphase which similar as for  Kay style relaxation
                            % 1.85 A corresponds to 13.7 / 24.7 inphase/antiphase which is similar to as reported by Jiang et al.
phiN    = 20*pi/180;        % angle dipolar and CSA frames for 15N [rad]
phiH    = 0;                % angle dipolar and CSA frames for 1H [rad]
mu0f    = 1e-7;             % permeability of vacuum factor mu0/(4pi)
                            % Palmer: mu0 = 4*pi*1e-7
ti      = 1/(1/tc + 1/te);  % internal motion effective correlation time Allard Eq. 16

% derive offsets
wN = 1e-6*gN*B0*wNppm;     % offset N in rad s-1
wH = 1e-6*gH*B0*wHppm;     % offset H in rad s-1

% derive input for LVM

offset = offResonance*[wH wN];
gB1_H  = 2*pi*1./(4*pw_dec)*1e6;
gB1_N  = 2*pi*1/(4*pw_90)*1e6;
gB1_cp = 2*pi*1/(4*pw_cpmg)*1e6;


% derived dispersion parameters
numPulses   = [2:2:2*time_T2*maxCPMG];        % minimum numpulses is 2, oner per CPMG half, up to 1000 Hz which 2*time_T2 in ms, in steps of 2
                                              % note that numPulses is split in two halves w/ P-element in middle
npoints     = length(numPulses);              % number of points in dispersion curve
nuCPMG      = 0.5*numPulses/time_T2;          % 2*tau between pulses so nuCPMG=1/4*tau  AND 2*k*2*tau=CT <=> tau = CT/(4k) = CT/(2N)
                                              % thus nuCPMG = 1/(4*CT/2N) = 0.5*N/CT

if feedback >= 1
    disp("")
    disp(" *** System parameters ***")
    printf("\tOffset 15N: %6.1f Hz\n", wNppm*gN*B0/(2*pi)*1e-6)
    printf("\tOffset  1H: %6.1f Hz\n", wHppm*gH*B0/(2*pi)*1e-6)
    printf("\tgB2     1H: %6.1f kHz\n", gB1_H/(2*pi)*1e-3)
end

% define start magnetization Ny (in-phase 15N as in CW experiment, ST-CW, ST-CPD)
StartOp     = zeros(16,1);
StartOp(7)  = 1;              % Nz=7 Ny=6 in Allard style matrix
StartOp(4)  = 1;              % Hz=4 in Allard style matrix
StartOp(1)  = 0.5;            % Em in Allard style matrix
