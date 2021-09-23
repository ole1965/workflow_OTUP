function stop = outfun_cluster(x,optimValues,state)
    stop = false;
    if exist('rfPulses\currentPulses_cluster\currentPulse.mat')
        load rfPulses\currentPulses_cluster\currentPulse.mat 
        zeit2 = datetime;
        iters = [iters;optimValues.iteration];
        zeit = [zeit;datetime];
        pulseSave = [pulseSave,x];
        fcnvalues = [fcnvalues,optimValues.fval];
        save rfPulses\currentPulses_cluster\currentPulse.mat iters zeit pulseSave fcnvalues zeit2
    else
        zeit2 = datetime;
        iters = [optimValues.iteration];
        zeit = [datetime];
        pulseSave = [x];
        fcnvalues = [optimValues.fval];
        save rfPulses\currentPulses_cluster\currentPulse.mat iters zeit pulseSave fcnvalues zeit2
    end
end