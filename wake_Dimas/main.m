clear all;
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'n';
eigspec = 'all';
N = 1000;
k = 3;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
tic;
[A, B] = makeAB(k, N, h(k), Fr);
if strcmpi(bal,'y')
    [o, an] = balanceAB(A, B, eigspec, alg);
else
    [o, an] = solveGEP(A, B, eigspec, alg);
end
if  strcmpi(eigspec,'all')
    c = o./k;
    o1 = k.*filt(c);
    plot(real(o),imag(o),'o');
else
    o1 = o;
end
toc;

%% Plot eigenvalue spectrum
figure;
scatter(real(o),imag(o),'ok');
hold on;
yline(0,'linewidth',1.5,'color','#898989');
scatter(real(o1),imag(o1),'b','filled');
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);
