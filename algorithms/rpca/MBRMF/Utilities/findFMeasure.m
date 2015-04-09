function [fmax S bestTH] = findFMeasure(E, S0)
    fmax = 0;
    Vth = linspace(min(E(:)), max(E(:)), 10000);
    for idx = 1:length(Vth)
        th = Vth(idx);
        Stmp = abs(E) > th;
        [pre rec] = countPR(S0, Stmp);
        ftmp = 2*rec*pre/(pre+rec);
        if ~isnan(ftmp) && ftmp > fmax
            fmax = ftmp;
            S = Stmp;
            bestTH = th;
        end
    end
end