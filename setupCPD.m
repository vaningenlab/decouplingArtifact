% setupCPD.m

cpdString = {"CW", "TWFO", "MLEV","WALTZ", "SPA"};

if cpd_type == 0        % CW
    cpd_phase = [ 0   0   0];                                        % phase of each element 0=x 1=y 2=-x 3=-y
    cpd_angle = [90 240  90];                                        % pulse flip angles of each element
elseif cpd_type == 3    % WALTZ
    cpd_phase = [ 0   2   0];                                        % phase of each element 0=x 1=y 2=-x 3=-y
    cpd_angle = [90 180 270];                                        % pulse flip angles of each element
elseif cpd_type == 2   % MLEV
    cpd_phase = [ 0   1   0];                                        % phase of each element 0=x 1=y 2=-x 3=-y
    cpd_angle = [90 180  90];                                        % pulse flip angles of each element
elseif cpd_type == 1   % 240
    cpd_phase = [ 0   1   0];                                        % phase of each element 0=x 1=y 2=-x 3=-y
    cpd_angle = [90 240  90];                                        % pulse flip angles of each element
elseif cpd_type == 4   % 27-126-27 Korolova et al 2013
    cpd_phase = [ 3   1  3 ];                                        % phase of each element 0=x 1=y 2=-x 3=-y
    cpd_angle = [27 126  27 ];                                        % pulse flip angles of each element
else
    disp("\nAborting, unkown type of CPD block...\n")
    exit
end

% check CPD
if length(cpd_phase) > 3 || length(cpd_angle) > 3
    disp("\nAborting, more than 3 elements in CPD is not supported...\n")
    exit
end

if length(cpd_phase) != length(cpd_angle)
    disp("\nAborting, CPD phase and angle arrays have different sizes...\n")
    exit
end

% derive CPD element, e.g. [0 1 0 0] for element 1 with [wHx wHy wNx wNy] notation
for c=1:length(cpd_phase)
    if cpd_phase(c) == 0 cpdPhase(c,:) = [1 0 0 0]; end
    if cpd_phase(c) == 1 cpdPhase(c,:) = [0 1 0 0]; end
    if cpd_phase(c) == 2 cpdPhase(c,:) = [-1 0 0 0]; end
    if cpd_phase(c) == 3 cpdPhase(c,:) = [0 -1 0 0]; end
end
