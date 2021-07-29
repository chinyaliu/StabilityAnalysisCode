clear all;
if ~contains(path,'code_Dimas;')
    addpath('code_Dimas');
end    
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'y';
eigspec = 'max';
N = 1500;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
klist = linspace(0.01,4,400);
if  strcmpi(eigspec,'all')
    olist = cell(length(klist),1);
    o_dis = cell(length(klist),1);
else
    olist = nan(1,length(klist));
end
tic;
for i = 1:length(klist)
    k = klist(i);
    fprintf('k = %.2f\n', k);
    [A, B] = makeAB(k, N, h(k), Fr);
    if strcmpi(bal,'y')
        o = balanceAB(A, B, eigspec, alg);
    else
        o = solveGEP(A, B, eigspec, alg);
    end
    if  strcmpi(eigspec,'all')
        o1 = k.*filt(o./k);
        olist{i} = o;
        o_dis{i} = o1;
    else
        if ~isempty(o)
            olist(i) = o(1);
        end
    end
end
toc;

%% Plot eigenvalue spectrum
figure;
if  strcmpi(eigspec,'all')
    hold on;
    for i = 1:length(olist)
        plot(real(o_dis{i}),imag(o_dis{i}),'ko');
    end
    hold off;
    xlim([0 1.4]);
    ylim([-0.6 0.1]);
else
    plot(real(olist),imag(olist),'ko');
end

%% Plot growth rate
figure;
plot(real(klist),imag(olist),'-k','linewidth',2);
xlabel('$k$');
ylabel('$\omega_i$');
grid on;