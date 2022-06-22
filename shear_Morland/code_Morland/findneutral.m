function [nmin, nmax] = findneutral(ud,dlt)
% Find neutral stability modes for EXPONENTIAL profile
    polylam = [4*ud^2+pi*dlt 0 -4*pi*dlt*ud^4+4*ud^2+2*pi*dlt 0 pi*dlt];
    lamall = roots(polylam);
    lam = real(lamall(imag(lamall)<eps & real(lamall)>eps));
    if length(lam) == 2
        nmin = min(lam);
        nmax = max(lam);
    elseif length(lam) == 1
        nmin = lam;
        nmax = NaN;
    else
        nmin = NaN;
        nmax = NaN;
    end
end