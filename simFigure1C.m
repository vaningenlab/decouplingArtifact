% simAll.m

% this is a simulation in a simple 2-spin space exactly as described by Allard, not usable for exchange!
% this way was also use in Jiang et al.

ini                         % initialize

B0array = [ 500 600 700 750 800 850 900 950 1200];

for bb=1:length(B0array)
    defineParsFig1C             % define basic parameters <++ ADJUST TO YOUR LIKINGS ++>
    buildRelaxationMatrix       % derive relation rates and build matrix
    for oo=1:length(wHppm)
        if feedback >= 1
            disp("")
            printf("*** B0: %6d MHz // 1H offset %6.1f ppm ***\n", B0array(bb), wHppm(oo))
        end
        % adjustment for calculation of Figure 1
        offset = offResonance*[wH(oo) wN];
        CPMG_CW                      % Flemming as well
        CPMG_ST_CW                   % now do simulation of Jiang's 15N CPMG with a single strain CW at fixed strength and Zuiderweg phase cycle
        R2c_ST    = R2_ST + ((laN-rhoN)*pw_cpmg*1*1e-6)./(1./nuCPMG);
        R2inf_ST     = sum(R2c_ST(kk/2:kk))/(0.5*kk+1);
        R2inf_CW     = sum(R2_CW(kk/2:kk))/(0.5*kk+1);
        maxDev_ST(oo,bb) = max(abs(R2c_ST - R2inf_ST));
        maxDev_CW(oo,bb) = max(abs(R2_CW - R2inf_CW));
    end
end

figure(1)
hold off
%plot(B0array,maxDev_CW(1,:),'bo-;1.5 ppm CW;')
plot(B0array,maxDev_CW(1,:),'bo-;1.0 ppm CW;')
hold on
%plot(B0array,maxDev_CW(2,:),'ro-;3.0 ppm CW;')
plot(B0array,maxDev_CW(2,:),'ro-;2.0 ppm CW;')

%plot(B0array,maxDev_ST(1,:),'b*-;1.5 ppm ST;')
plot(B0array,maxDev_ST(1,:),'b*-;1.0 ppm ST;')
%plot(B0array,maxDev_ST(2,:),'r*-;3.0 ppm ST;')
plot(B0array,maxDev_ST(2,:),'r*-;2.0 ppm ST;')

xlabel("B0 (MHz)")
ylabel("maximum deviation in R2,eff (s-1)")
plot([475 1225], [0.3 0.3], 'm')
grid on
axis([475 1225 -0.1 2.0])


if feedback >= 1
    disp("")
    toc                         % report timing
    disp("")
end
