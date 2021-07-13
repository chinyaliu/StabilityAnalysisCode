clear all;
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'n';
eigspec = 'all';
N = 1000;
% k = 0.6;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
klist = linspace(0.01,4,100);
olist = nan(1,length(klist));
tic;
for i = 1:length(klist)
    k = klist(i);
    fprintf('k = %.2f\n', k);
    [A, B] = makeAB(k, N, h(k), Fr);
    if strcmpi(bal,'y')
        [o, an] = balanceABs(A, B, eigspec, alg);
    else
        [o, an] = solveGEPs(A, B, eigspec, alg);
    end
    if  strcmpi(eigspec,'all')
        c = o./k;
        o1 = k.*filt(c);
    else
        o1 = o;
    end
    if ~isempty(o1)
        olist(i) = o1(1);
    end
end
toc;

%% Plot eigenvalue spectrum
figure;
plot(real(olist),imag(olist),'o');

%% Plot growth rate
figure;
plot(klist,imag(olist),'o');
