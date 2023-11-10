function compassDir = vectorInCompassDir(thisdx, thisdy)
    
    thislength = sqrt(thisdx^2 + thisdy^2);
    ThisSin = thisdy/thislength;
    ThisCos = thisdx/thislength;
    if ((ThisCos >= 0)&&(ThisSin == 0)) 
        ThisDir = 90;
    elseif ((ThisCos < 0)&&(ThisSin == 0)) 
        ThisDir = 270;    
    elseif (ThisCos < 0) 
        ThisDir = 180 + 90 - atand(ThisSin/ThisCos); 
    elseif (ThisCos > 0) 
        ThisDir = 90 - atand(ThisSin/ThisCos);
    elseif ThisCos == 0 && ThisSin > 0
        ThisDir = 0;
    elseif ThisCos == 0 && ThisSin < 0
        ThisDir = 180;
    end
    compassDir = mod(ThisDir + 360, 360);
end