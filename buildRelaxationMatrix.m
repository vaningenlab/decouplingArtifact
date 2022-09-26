% calculate relaxation rates

% HvI 2016 changed from Palmer to Allard style

    wHr = gH*B0;                                                                   % 1. evaluate resonance frequencies (ignoring chemical shift)
    wNr = gN*B0;
    J0    = 2/5*( S2*tc/1 + (1-S2)*ti/1);                                          % 2. define spectral densities
    JwH   = 2/5*( S2*tc/(1+(wHr*tc)^2) + (1-S2)*ti/(1+(wHr*ti)^2));
    JwN   = 2/5*( S2*tc/(1+(wNr*tc)^2) + (1-S2)*ti/(1+(wNr*ti)^2));
    JwHpN = 2/5*( S2*tc/(1+((wHr+wNr)*tc)^2) + (1-S2)*ti/(1+((wHr+wNr)*ti)^2));
    JwHmN = 2/5*( S2*tc/(1+((wHr-wNr)*tc)^2) + (1-S2)*ti/(1+((wHr-wNr)*ti)^2));
    JwHpH = 2/5*( S2*tc/(1+((wHr+wHr)*tc)^2) + (1-S2)*ti/(1+((wHr+wHr)*ti)^2));
                                                                                   % 3. define pre-factors
    dHN = 3*mu0f*hbar*gN*gH*rNH^(-3);   % DD interaction constant: Helgstrand/Allard has factor 3! see below Eq. 14
    cN = -CSAN*gN*B0;                   % CSA interaction constant 15N; Eq. 15
    cH = -CSAH*gH*B0;                   % CSA interaction constant 1H; not used by Allard, only by Helgstrand
                                                                                    % 4. relaxation due to external proton close to amide proton
    dHH   = 3*mu0f*hbar*gH*gH*rHH^(-3); % DD interaction constant: Helgstrand/Allard has factor 3! see below Eq. 14
    laHH  = 1/36*dHH^2*( 2*J0 + 3/2*JwH + 1/2*J0 + 3*JwH + 3*JwHpH);  % this gives reasonable values comped to Kay style for proton at 2.25 A
    rhoHH = 1/36*dHH^2*( 3*JwH + 1*J0 + 6*JwHpH);
                                                                                   % 5.longitudinal rates
    rhoN   = 1/36*dHN^2*(3*JwN + 1*JwHmN + 6*JwHpN) + 1/3*cN^2*JwN;                % Nz auto Eq. 25
    rhoH   = 1/36*dHN^2*(3*JwH + 1*JwHmN + 6*JwHpN) + 1/3*cH^2*JwH + rhoHH;        % Hz auto Eq. 26
    rhoHN  = 1/36*dHN^2*(3*JwN + 3*JwH) + 1/3*cN^2*JwN + 1/3*cH^2*JwH + rhoHH;     % HzNz auto Eq. 30 / Helgstrand 31
    sigma  = 1/36*dHN^2*(6*JwHpN - 1*JwHmN);                                       % cross-relaxation Eq. 31
    deltaN = 1/3*dHN*cN*1/2*(3*(cos(phiN))^2-1)*JwN;                               % cross-correlation Eq. 33
    deltaH = 1/3*dHN*cH*1/2*(3*(cos(phiH))^2-1)*JwH;                               % cross-correlation Helgstrand Eq. 35
                                                                                   % 6. transverse rates
    laN = 1/36*dHN^2*( 2*J0 + 3/2*JwN + 1/2*JwHmN + 3*JwH + 3*JwHpN) ...
            + 1/3*cN^2*(2/3*J0+1/2*JwN);                                      % Nxy auto Eq.23
    laH = 1/36*dHN^2*( 2*J0 + 3*JwN + 1/2*JwHmN + 3/2*JwH + 3*JwHpN) ...
            + 1/3*cH^2*(2/3*J0+1/2*JwH) + laHH;                               % Hxy auto Eq.24
    rhoaN = 1/36*dHN^2*(2*J0 + 3/2*JwN + 1/2*JwHmN + 3*JwHpN) ...
            + 1/3*cH^2*JwH + 1/3*cN^2*(2/3*J0 + 1/2*JwN) + rhoHH;             % HzNxy auto Eq. 27 / Helgstrand Eq. 28
    rhoaH = 1/36*dHN^2*(2*J0 + 3/2*JwH + 1/2*JwHmN + 3*JwHpN) ...
            + 1/3*cN^2*JwN + 1/3*cH^2*(2/3*J0 + 1/2*JwH) + laHH;              % HxyNz auto Eq. 28 / Helgstrand Eq. 29
    etaN = 1/3*dHN*cN*1/2*(3*(cos(phiN))^2-1)*(2/3*J0+1/2*JwN);                    % cross-correlation Eq. 34 //// CORRECTED TYPO! factor 1/2 was missing /////
    etaH = 1/3*dHN*cH*1/2*(3*(cos(phiH))^2-1)*(2/3*J0+1/2*JwH);                    % cross-correlation Helgstrand Eq. 37 //// CORRECTED TYPO! factor 1/2 was missing /////
                                                                                   % 7. multiple quantum rates
    laMQ = 1/36*dHN^2*(3/2*JwN + 1/2*JwHmN + 3/2*JwH + 3*JwHpN) ...
            + 1/3*cN^2*(2/3*J0 + 1/2*JwN) + 1/3*cH^2*(2/3*J0 + 1/2*JwH) + laHH;  % Eq. 29 / Helgstrand Eq. 30
    muMQ = 1/36*dHN^2*(3*JwHpN - 1/2*JwHmN);                                         % Eq. 32.
                                                                                % 7. define factors for unity matrix
    MH0 = 1;                        % Allard: MI0 and MS0 are the equilibrium magnetizations of spin I and S, respectively
    MN0 = 1*gN/gH;                  % Helgstrand (w/ exchange Eq. 18 and 19, scaled for fractional population of state)
    thH = rhoH*MH0 + sigma*MN0;     % Eq. 38 Helgstrand
    thN = rhoN*MN0 + sigma*MH0;     % Eq. 39 Helgstrand
    thHN = deltaH*MH0 + deltaN*MN0; % Eq. 40 Helgstrand
    
    relaxRates = includeRelaxation*[laH etaH thH rhoH sigma deltaH laN etaN thN rhoN deltaN rhoaH rhoaN laMQ muMQ thHN rhoHN];
    
if feedback >= 1
    disp("")
    disp(" *** Calculated relaxation rates ***")
    printf("\tR2 Nxy   is %3.3f\n", laN)
    printf("\tR2 NxyHz is %3.3f\n", rhoaN)
    printf("\tR1 Nz    is %3.3f\n", rhoN)
    printf("\tR1 NzHz  is %3.3f\n", rhoHN)
    printf("\tCCR Nxy  is %3.3f\n", etaN)
    printf("\tCCR Nz   is %3.3f\n", deltaN)
    printf("\tHz<->Nz  is %3.3f\n", sigma)
    printf("\tR2 Hxy   is %3.3f\n", laH)
    printf("\tR1 Hz    is %3.3f\n", rhoH)
end
