% simCPD.m

% this is to compare 240/MLEV/WALTZ performance using the matching setup
% and using the actual implementation of start point of CW/CPD so that there is no need for special matching of first spin-echo

ini                         % initialize
defineParsCPD               % define basic parameters <++ ADJUST TO YOUR LIKINGS ++>
buildRelaxationMatrix       % derive relation rates and build matrix

CPMG_CW                     % simulation of Flemming's 15N CPMG with matched 1H decoupling strength and without Zuiderweg phase cycle -- reference data
CPMG_ST_CW                  % now do simulation of Jiang's 15N CPMG with a single strain CW at fixed strength and Zuiderweg phase cycle -- reference data

for cc = 0:4                %   
    cpd_type  = cc;         % 0 = CW; 3 = WALTZ; 2= MLEV; 1=90-240-90; 4= 26-y 127+y 26-y = SPA
    setupCPD;
    CPMG_ST_CPD             % now do simulation of our 15N CPMG with a single strain CPD-train at matched strength and Zuiderweg phase cycle
    R2_CPD(cc+1,:)   = R2_CPDm;
    R2c_CPD(cc+1,:)  = R2_CPDm + ((laN-rhoN)*pw_cpmg*1*1e-6)./(1./nuCPMG);
    R2a_CPD(cc+1)    = sum(R2c_CPD(cc+1,:))/length(R2c_CPD(cc+1,:));
    RMSD_CPD(cc+1)   = sqrt(sum((R2c_CPD(cc+1,:) - R2a_CPD(cc+1)).^2)/length(R2c_CPD(cc+1,:)));
    maxDev_CPD(cc+1) = max(abs(R2c_CPD(cc+1,:) - R2a_CPD(cc+1)));
end

makePlotCPD                % plot results, only R1-corrected data for ST exps.

if feedback >= 1
    disp("")
    toc                    % report timing
    disp("")
end
