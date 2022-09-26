% mymod.m

function t=mymod(arg1, arg2)

    t = mod(arg1,arg2);
    if t==0
        t = arg2;
    end
    
endfunction

