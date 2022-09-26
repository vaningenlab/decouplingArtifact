% makePlot.m
oo=1;
plotColors = 'kgrcm'; 
if plot_hold == 1
    if oo==1
        hold off
        colorIdx = 1;
        plotNum = 1;
    else
        hold on
    end
    if mod(plotNum, 4) == 0
        colorIdx = 4;
    else
        colorIdx = mod(plotNum, 4);
    end
else
    hold off
    plotNum = 1;
    plotMax = 0;
    plotMin = 0;
    colorIdx = 1;
end

R2_STc   = R2_ST + ((laN-rhoN)*pw_cpmg*1*1e-6)./(1./nuCPMG);
figure(1)
hold off
plot(nuCPMG, R2_CW,'bo;CW;')            % Flemming
hold on
plot(nuCPMG, R2_STc,'r*-.;ST;')         % Jiang
plot(nuCPMG, R2_CWf,'m+-;CWf;')         % Flemming fixed

% add theoretical curve from Jiang
Jr     = JNH*wH / sqrt( gB1_H^2 + wH^2);
tau_cp = 1./(4*nuCPMG);
% this is now corrected -- factor 2
R2_cal = R2_CW(npoints) + 0.5*(rhoaN-laN)*(1 - sinc(2*Jr*tau_cp) );
plot(nuCPMG, R2_cal, 'ks--;eq. 1;' )        % equation 1 based 

title("N15 relaxation dispersion curves")
xlabel("vCPMG (Hz)")
ylabel("R2,eff (s-1)")
grid on
axis([0 1010 13.5 16.5])

if feedback >= 1
    disp("")
    disp(" *** size slow pulsing artifact ***")
    printf("\tCW   : %6.1f s-1\n", R2_CW(1)-R2_CW(npoints))
    printf("\tCWf  : %6.1f s-1\n", R2_CWf(1)-R2_CWf(npoints))
    printf("\tST-CW: %6.1f s-1\n", R2_STc(1)-R2_STc(npoints))
    printf("\tEq. 1: %6.1f s-1\n", R2_cal(1)-R2_cal(npoints))
end


