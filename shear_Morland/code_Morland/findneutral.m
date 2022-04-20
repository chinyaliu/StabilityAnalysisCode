function out = findneutral(ud,dlt)
% Find the neutral stability curve of exponential single shear profile.
% Parameters are given in dimensionless form
% Velocity profile U(z) = ud*exp(2z/dlt)
% Neutral curve [4*ud^2+pi*dlt](lam^4) + [-4pi*dlt*ud^4+4*ud^2+2pi*dlt](lam^2) + pi*dlt = 0
%
fx = @(x,y) [4*x^2+pi*y 0 -4*pi*y*x^4+4*x^2+2*pi*y 0 pi*y];
out = cell(1,length(dlt));
for i = 1:length(dlt)
    polylam = fx(ud,dlt(i));
    lamm = roots(polylam);
    out{i} = real(lamm(imag(lamm)<eps & real(lamm)>-eps));
end
end