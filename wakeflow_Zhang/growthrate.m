close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4", "Ray_match"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];%, "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 400;
k = linspace(0.01,4,400);
Re = inf;
Fr2 = 2.25;
inflec_pt = -0.74708299;
cutz(1) = -inflec_pt;
h1 = @(k) 2*pi/k;
h2 = @(k) 6.5;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
tic;
p1 = wZhang_solver(N,1,1,Re,Fr2,method);
for i = 1:length(k)
    if (k(i) < pi/3)
        h = h1;
    else
        h = h2;
    end
    p1.k = k(i); p1.h = h(k(i));
    [o(i), ~, ~, ~, ~] = p1.solver(cutz(i), 'y', alg, do_balancing);
    z(:,i) = p1.z; phi{i} = p1.phi; z_c(i) = p1.zc; cutz(i+1)=p1.zL;
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
toc;
cutz = cutz(2:end);
% save('modeshape','phi','z','N','k','cutz','o','z_c');
%% Plot growth rate v.s. k
o1 = imag(o);
o1r = real(o);
fig1 = figure('position',[50,0,1000,720]);
plot(k,o1,'linewidth',2);
ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;