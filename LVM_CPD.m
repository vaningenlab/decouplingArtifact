% LVM_CPD.m
% function to propagate magnetization in presence of CPD block on 1H
% version to be used when matching the power level
% version adapted to do totalDuration shorter than tp_CPD

% needs to be able to call LVM
%       thus needs to know offset, JNH, gB1_H/gB1_N, relaxRates !
%       these are passed as input inside a container

function cv_  = LVM_CPD( inputMagnetisation , totalDuration, remainingBlockTime, nitrogenPulse, lvPars)

    global remTime                           % remTime will be modified and thus need be inside function-loop
    
    global CPD1 CPD2 CPD3 CPD1_Nyp CPD2_Nyp CPD3_Nyp CPD1_Nxp CPD2_Nxp CPD3_Nxp CPD1_Nxm CPD2_Nxm CPD3_Nxm
    global CPDt CPDt_Nyp CPDt_Nxp CPDt_Nxm
    global V1_0 V2_0 V3_0 V1_Nyp V2_Nyp V3_Nyp V1_Nxp V2_Nxp V3_Nxp V1_Nxm V2_Nxm V3_Nxm
    global D1_0 D2_0 D3_0 D1_Nyp D2_Nyp D3_Nyp D1_Nxp D2_Nxp D3_Nxp D1_Nxm D2_Nxm D3_Nxm
    global tp_1 tp_2 tp_3 tp_CPD

    % shorthands
    inpMag  = inputMagnetisation;
    totDur  = totalDuration;
    remTime = remainingBlockTime;
    gB1_N_  = nitrogenPulse;     % contains strength and axis of 15N rotation
    
    % extract constants (within one simulation)
    offset_ = lvPars{1};
    piJ_    = lvPars{2};
    relaxR_ = lvPars{3};
    cpdAngle= lvPars{4};
    
    cv_    = inpMag;            % to keep code similar to main body
    elapTime = 0;               % to keep track of how much time has actually elapsed
    
    % find out which pre-defined matrix-exponential / V,D-pair to use for full-element/partial element propagators
    if gB1_N_ == 0
        LV_CPD1  = CPD1;
        LV_CPD2  = CPD2;
        LV_CPD3  = CPD3;
        LV_CPDf  = CPDt;
        V1    = V1_0;
        D1    = D1_0;
        V2    = V2_0;
        D2    = D2_0;
        V3    = V3_0;
        D3    = D3_0;
    elseif gB1_N_(4) > 0
        LV_CPD1  = CPD1_Nyp;
        LV_CPD2  = CPD2_Nyp;
        LV_CPD3  = CPD3_Nyp;
        LV_CPDf  = CPDt_Nyp;
        V1       = V1_Nyp;
        D1       = D1_Nyp;
        V2       = V2_Nyp;
        D2       = D2_Nyp;
        V3       = V3_Nyp;
        D3       = D3_Nyp;
    elseif gB1_N_(3) > 0
        LV_CPD1  = CPD1_Nxp;
        LV_CPD2  = CPD2_Nxp;
        LV_CPD3  = CPD3_Nxp;
        LV_CPDf  = CPDt_Nxp;
        V1       = V1_Nxp;
        D1       = D1_Nxp;
        V2       = V2_Nxp;
        D2       = D2_Nxp;
        V3       = V3_Nxp;
        D3       = D3_Nxp;
    elseif gB1_N_(3) < 0
        LV_CPD1  = CPD1_Nxm;
        LV_CPD2  = CPD2_Nxm;
        LV_CPD3  = CPD3_Nxm;
        LV_CPDf  = CPDt_Nxm;
        V1       = V1_Nxm;
        D1       = D1_Nxm;
        V2       = V2_Nxm;
        D2       = D2_Nxm;
        V3       = V3_Nxm;
        D3       = D3_Nxm;
    else
        printf("\n ERROR input of nitrogenPulse!\n")
    end
    
    % apply the remainder of the block that was not finished before
    % remainingBlockTime needs always less than a complete blockTime (90+240+90) by definition
    % so we can judge where to pick up the CPD within the block based on the length of remainingBlockTime
    
    if remTime > tp_CPD
        printf("\n**TIMING ERROR** remaining time is %.2f (us), more than one CPD-block %.2f (us)!!\n", 1e6*remTime, 1e6*tp_CPD)
        printf("**-- ABORTED --**")
        exit
    end
    
    if remTime > (tp_2 + tp_3)
        pulseTime = min(remTime - (tp_2 + tp_3), totDur);
        LVexp     = real(V1*diag(exp(diag(D1)*pulseTime))/V1); % works maybe need to take real part
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;       % new elapsed time
        pTimeRem  = max(totDur-elapTime, 0);    % total time remaining to execute this function
        pulseTime = min(tp_2, pTimeRem);        % new pulse time for next element
        LVexp     = real(V2*diag(exp(diag(D2)*pulseTime))/V2);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;       % new elapsed time, note pulseTime can be zero
        pTimeRem  = max(totDur-elapTime, 0);    % total time remaining to execute this function
        pulseTime = min(tp_3, pTimeRem);        % new pulse time for next element
        LVexp     = real(V3*diag(exp(diag(D3)*pulseTime))/V3);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;        % new elapsed time, note pulseTime can be zero
    elseif remTime > tp_3
        % first apply remainder of left-over central element (this could be complete)
        pulseTime = min(remTime - (tp_3), totDur);
        LVexp     = real(V2*diag(exp(diag(D2)*pulseTime))/V2);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;       % new elapsed time
        pTimeRem  = max(totDur-elapTime, 0);    % total time remaining to execute this function
        pulseTime = min(tp_3, pTimeRem);        % new pulse time for next element
        LVexp     = real(V3*diag(exp(diag(D3)*pulseTime))/V3);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;        % new elapsed time, note pulseTime can be zero
    else
        % first apply remainder of left-over last element (this could be complete)
        pulseTime = min(remTime, totDur);
        LVexp     = real(V3*diag(exp(diag(D3)*pulseTime))/V3);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;
    end
    % check timing
    if abs(elapTime-remTime) > 0.001e-9 % 1 ns accuracy
        if totDur > remTime
            printf("\n**TIMING ERROR** remainingBlockTime was %g (s), elapsed time is %g (s)!!\n", remTime, elapTime)
        else
            if abs(elapTime-totDur) > 0.001e-9
                printf("\n**TIMING ERROR** total duration was %g (s), elapsed time is %g (s)!!\n", totDur, elapTime)
            end
        end
    end
    
    % figure out how many full CPD blocks we can do in the remainder of available time
    loopTime = max(totDur - elapTime, 0);
    if loopTime > tp_CPD            % else we can skip this and only do the last part
        numBlocks = floor(loopTime/tp_CPD);
        for bb=1:numBlocks
            cv_ = LV_CPDf*cv_;
        end
        elapTime = elapTime + numBlocks*tp_CPD;
    end
    
    % figure out how much of a block is left to do
    timeLeft = totDur - elapTime;
    if timeLeft >= tp_CPD
        printf("\n**TIMING ERROR** timeLeft is %g (s), bigger than/equal to time for a CDP block, which is %g (s)!!\n", timeLeft, tp_CPD)
    end
    if timeLeft < 0
        printf("\n**TIMING ERROR** timeLeft is %g (us), negative!!\n", 1e6*timeLeft)
    end
    if timeLeft > (tp_2 + tp_1)
        % first apply element 1 and 2 then remainder of final element
        cv_       = LV_CPD1*cv_;
        cv_       = LV_CPD2*cv_;
        pulseTime = timeLeft - (tp_2 + tp_1);
        LVexp     = real(V3*diag(exp(diag(D3)*pulseTime))/V3);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + tp_1 + tp_2 + pulseTime;
        remTime   = tp_CPD   - tp_1 - tp_2 - pulseTime;
    elseif timeLeft > tp_1
        % first apply element 1 then remainder of central element (this could be complete)
        cv_       = LV_CPD1*cv_;
        pulseTime = timeLeft - (tp_1);
        LVexp     = real(V2*diag(exp(diag(D2)*pulseTime))/V2);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + tp_1 + pulseTime;
        remTime   = tp_CPD   - tp_1 - pulseTime;
    elseif timeLeft > 0
        % first apply remainder of initial element (this could be complete)
        pulseTime = timeLeft;
        LVexp     = real(V1*diag(exp(diag(D1)*pulseTime))/V1);
        cv_       = LVexp*cv_;
        elapTime  = elapTime + pulseTime;
        remTime   = tp_CPD   - pulseTime;
    end
    rt_ = remTime;
    %printf("OLD remTime is %.3f us, totalDuration is %.3f us, blockLength is %.3f us, blocks is %.3f\n", remainingBlockTime*1e6, totalDuration*1e6, tp_CPD*1e6, totalDuration/tp_CPD)
    %printf("NEW remTime is %.3f us\n", rt_*1e6)
    % check timing
    if abs(elapTime-totDur) > 0.001e-9      % 1ns accuracy
        printf("\n**TIMING ERROR** totalDuration was %f (s), elapsed time is %f (s)!!\n", totDur, elapTime)
    end
    % check time remaining
    if remTime > tp_CPD
        printf("\n**TIMING ERROR** remaining time is %.2f (us), more than one CPD-block %.2f (us)!!\n", 1e6*remTime, 1e6*tp_CPD)
        printf("\n\n timeLeft was  %.2f (us), pulseTime %.2f, tp_1 / tp_2 is %.2f\n", 1e6*timeLeft, 1e6*pulseTime, 1e6*tp_1,1e6*tp_2)
    end

endfunction
