% LVM.m

function m = LVM(offset, piJ, pulse, relaxRates)

    % process input
    wH  = offset(1);
    wN  = offset(2);
    wHx = pulse(1);
    wHy = pulse(2);
    wNx = pulse(3);
    wNy = pulse(4);
    laH     = relaxRates(1);
    etaH    = relaxRates(2);
    thH     = relaxRates(3);
    rhoH    = relaxRates(4);
    sigma   = relaxRates(5);
    deltaH  = relaxRates(6);
    laN     = relaxRates(7);
    etaN    = relaxRates(8);
    thN     = relaxRates(9);
    rhoN    = relaxRates(10);
    deltaN  = relaxRates(11);
    rhoaH   = relaxRates(12);
    rhoaN   = relaxRates(13);
    laMQ    = relaxRates(14);
    muMQ    = relaxRates(15);
    thHN    = relaxRates(16);
    rhoHN   = relaxRates(17);

    %exact formulation of Allard
    %        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16 
    %        E/2  Ix   Iy   Iz   Sx   Sy   Sz   IxSz IySz IzSx IzSy IxSx IxSy IySx IySy IzSz
    m  = -1*[  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    ;
            0    laH  wH   -wHy 0    0    0    etaH piJ  0    0    0    0    0    0    0    ;
            0    -wH  laH  wHx  0    0    0    -piJ etaH 0    0    0    0    0    0    0    ;
           -thH  wHy  -wHx rhoH 0    0  sigma  0    0    0    0    0    0    0    0  deltaH ;
            0    0    0    0    laN  wN   -wNy 0    0    etaN piJ  0    0    0    0    0    ;
            0    0    0    0    -wN  laN  wNx  0    0    -piJ etaN 0    0    0    0    0    ;
           -thN  0    0  sigma  wNy  -wNx rhoN 0    0    0    0    0    0    0    0  deltaN ;
            0    etaH piJ  0    0    0    0  rhoaH  wH   0    0    wNy  -wNx 0    0    -wHy ;
            0    -piJ etaH 0    0    0    0  -wH  rhoaH  0    0    0    0    wNy  -wNx wHx  ;
            0    0    0    0    etaN piJ  0    0    0  rhoaN  wN   wHy  0    -wHx 0    -wNy ;
            0    0    0    0    -piJ etaN 0    0    0  -wN  rhoaN  0    wHy  0    -wHx wNx  ;
            0    0    0    0    0    0    0    -wNy 0  -wHy   0    laMQ wN   wH  -muMQ 0    ;
            0    0    0    0    0    0    0    wNx  0    0    -wHy -wN  laMQ muMQ wH   0    ;
            0    0    0    0    0    0    0    0    -wNy wHx  0    -wH  muMQ laMQ wN   0    ;
            0    0    0    0    0    0    0    0    wNx  0    wHx -muMQ -wH  -wN  laMQ 0    ;
          -thHN  0    0  deltaH 0    0  deltaN wHy  -wHx wNy  -wNx 0    0    0    0    rhoHN ];
        

endfunction
