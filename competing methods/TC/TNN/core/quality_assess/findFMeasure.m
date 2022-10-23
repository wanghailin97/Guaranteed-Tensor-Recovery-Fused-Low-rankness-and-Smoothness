function [fmax, S, bestTH] = findFMeasure(E, S0)
    fmax = 0;
    Vth = linspace(min(abs(E(:))), max(abs(E(:))), 5000);
    for idx = 1:length(Vth)
        th = Vth(idx);
        Stmp = abs(E) > th;
        [pre, rec] = countPR(S0, Stmp);
        ftmp = 2*rec*pre/(pre+rec);
        if ~isnan(ftmp) && ftmp > fmax
            fmax = ftmp;
            S = Stmp;
            bestTH = th;
        end
    end
    
    if fmax == 0
        S = Stmp;
        bestTH = 0;
    end
end

function         [pre, rec] = countPR(S0, Stmp)
CS = (S0.*Stmp);
pre = sum(CS(:))/sum(Stmp(:));
rec = sum(CS(:))/sum(S0(:));
end