% simFigure1A.m

% this is a simulation in a simple 2-spin space exactly as described by Allard, not usable for exchange!
% this way was also used in Jiang et al.

ini                         % initialize
defineParsFig1A             % define basic parameters <++ ADJUST TO YOUR LIKINGS ++>
buildRelaxationMatrix       % derive relation rates and build matrix

CPMG_CW                     % now do simulation of Flemming's 15N CPMG with CW at matched strength
CPMG_CW_fixed               % same but fixed CW strength
CPMG_ST_CW                  % now do simulation of Jiang's 15N CPMG with single train CW at fixed strength and XY16 phase cycle
makePlotFig1A               % plot results


if feedback >= 1
    disp("")
    toc                         % report timing
    disp("")
end
