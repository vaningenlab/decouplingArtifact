% simAll.m

% this is a simulation in a simple 2-spin space exactly as described by Allard, not usable for exchange!
% this way was also use in Jiang et al.

ini                         % initialize
defineParsFig1B             % define basic parameters <++ ADJUST TO YOUR LIKINGS ++>

for ff=1:length(tca)
    tc = tca(ff);
    ti      = 1/(1/tc + 1/te);  % internal motion effective correlation time Allard Eq. 16
    buildRelaxationMatrix       % derive relation rates and build matrix
    for oo=1:length(wHppm)
        if feedback >= 1
            disp("")
            printf("*** taucee: %6.1f ns // 1H offset %6.1f ppm ***\n", tc*1e9, wHppm(oo))
        end
        % adjustment for calculation of Figure 1
        offset = offResonance*[wH(oo) wN];
        CPMG_CW                      % Flemming as well
        CPMG_ST_CW                   % now do simulation of Jiang's 15N CPMG with a single strain CW at fixed strength and Zuiderweg phase cycle
        R2c_ST    = R2_ST + ((laN-rhoN)*pw_cpmg*1*1e-6)./(1./nuCPMG);
        R2inf_ST     = sum(R2c_ST(kk/2:kk))/(0.5*kk+1);
        R2inf_CW     = sum(R2_CW(kk/2:kk))/(0.5*kk+1);
        maxDev_ST(ff,oo) = max(abs(R2c_ST - R2inf_ST));
        maxDev_CW(ff,oo) = max(abs(R2_CW - R2inf_CW));
    end
end

figure(1)
hold off
plot(wHppm,maxDev_CW(1,:),'ko-;4 ns CW;')
hold on
plot(wHppm,maxDev_CW(2,:),'bo-;6.5ns CW;')
plot(wHppm,maxDev_CW(3,:),'ro-;9 ns CW;')
plot(wHppm,maxDev_ST(1,:),'k*-;4 ns ST;')
hold on
plot(wHppm,maxDev_ST(2,:),'b*-;6.5ns ST;')
plot(wHppm,maxDev_ST(3,:),'r*-;9 ns ST;')
xlabel("1H offset (ppm)")
ylabel("maximum deviation in R2,eff (s-1)")
plot([0 4], [0.3 0.3], 'm')
grid on
axis([0 4.1 -0.1 3])


if feedback >= 1
    disp("")
    toc                         % report timing
    disp("")
end
