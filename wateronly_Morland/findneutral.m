function [nmin, nmax] = findneutral(ud,dlt)
    polylam = [4*ud^2+pi*dlt 0 -4*pi*dlt*ud^4+4*ud^2+2*pi*dlt 0 pi*dlt];
    lamall = roots(polylam);
    lam = real(lamall(imag(lamall)<eps & real(lamall)>eps));
    if length(lam) == 2
        nmin = min(lam);
        nmax = max(lam);
    else
        nmin = NaN;
        nmax = NaN;
    end
end