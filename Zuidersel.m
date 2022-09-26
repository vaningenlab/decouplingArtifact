% mymod.m

function t=Zuidersel(arg1)

    p = mod(arg1,4);
    if p==1
        t = [0 0 0 1];      % y
    elseif p==2
        t = [0 0 0 1];      % y
    elseif p==3
        t = [0 0 1 0];      % x
    elseif p==0
        t = [0 0 -1 0];     % -x
    else
        printf("\n ERROR Zuiderweg selection!\n")
    end
    
endfunction

