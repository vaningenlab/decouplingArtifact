% makePlotCPD.m
labelShift=0.25;
labelWidth=0.1;

plotColors = 'kbrgmcy'; 
figure(1)
hold off
colorIdx = 1;

plot(nuCPMG, R2_CW,sprintf( '%ss-;CW;', plotColors(colorIdx) ))
colorIdx = colorIdx + 1;
hold on

plot(nuCPMG, R2_ST + ((laN-rhoN)*pw_cpmg*1*1e-6)./(1./nuCPMG),sprintf( '%so-;ST-CW;', plotColors(colorIdx) ))
colorIdx = colorIdx + 1;

text(200, min([R2_CW min(R2c_CPD)])-labelShift, "CPD     RMSD    max. deviation", 'fontweight', 'bold', 'fontsize', 11)
for cc=0:length(R2a_CPD)-1
    plot(nuCPMG, R2c_CPD(cc+1,:),sprintf( '%s*-;%s;', plotColors(colorIdx), cpdString{cc+1} ))
    text(200, min([R2_CW min(R2c_CPD)])-labelShift-(cc+1)*labelWidth, sprintf('%-6s %8.3f %6.2f', cpdString{cc+1}, RMSD_CPD(cc+1), maxDev_CPD(cc+1)),'fontsize', 11)
    colorIdx = colorIdx + 1;
end



if plot_jiang == 1
    axis([0 max(nuCPMG)+50 11.5 14.5]) 
else
    axis([0 max(nuCPMG)+50 round(min([R2_CW min(R2c_CPD)]))-1 round(max([R2_CW max(R2c_CPD)]))+1])
end

grid on
title("N15 relaxation dispersion curves")
xlabel("vCPMG (Hz)")
ylabel("R2,eff (s-1)")


