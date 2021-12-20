function [ENG,z] = energy_boomkamp(flow1,eval,evec)
[z, phi] = findmodeshape(flow1, evec);
[~, U] = getU(flow1);
k = flow1.k;
q = evec(end);

% KIN = DIS + REY + TEN + HYD + TAN
res = 0.5*k.*(imag(phi(:,2)).*real(phi(:,1))-real(phi(:,2)).*imag(phi(:,1)));
rey = res.*U(:,2);
kin = 0.5*imag(eval)*(imag(phi(:,2)).^2+real(phi(:,2)).^2+k^2*(real(phi(:,1)).^2+imag(phi(:,1)).^2));
dis = -0.5/flow1.Re*(4*k^2*(imag(phi(:,2)).^2+real(phi(:,2)).^2) ...
    +(real(phi(:,3))+k^2*real(phi(:,1))).^2 ...
    +(imag(phi(:,3))+k^2*imag(phi(:,1))).^2);
ten = 0; % ignore surface tension for wake problem
hyd = 0.5*k*(imag(q)*real(phi(1,1))-real(q)*imag(phi(1,1)))/flow1.Fr2;
tan = -0.5*U(1,3)*(real(q)*real(phi(1,2))+imag(q)*imag(phi(1,2)))/flow1.Re;

% Integrate from z = -h to z = 0
reyi = -trapz(z,rey);
kini = -trapz(z,kin);
disi = -trapz(z,dis);

% Scale by kinetic energy
[~,ind] = max(abs(kini));
scal = kini(ind);
REY = reyi/scal;
KIN = kini/scal;
DIS = disi/scal;
TEN = ten/scal;
HYD = hyd/scal;
TAN = tan/scal;
res = res/scal;
rey = rey/scal;
kin = kin/scal;
dis = dis/scal;

% Compute the relative error
tot = REY + DIS + TEN + HYD + TAN;
err = abs(tot-KIN);

ENG = struct('KIN',KIN,'DIS',DIS,'REY',REY,'TEN',TEN,'HYD',HYD,'TAN',TAN,...
    'total',tot,'error',err,'res',res,'rey',rey,'kin',kin,'dis',dis);
end