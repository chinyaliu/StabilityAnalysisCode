close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,2,3]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 20:20:800;
% k = 0.2;
% Re = 1e+9;
k = 1;
Re = 1e4;
%% Run solver with different grid numbers
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
for i = 1:length(N)
   [o, ~, cA(i), errG(i), ~] = poiseuille_solver(N(i),k,Re,method,alg,do_balancing);
   oi(i) = imag(o);
   fprintf("N = %3d, omega_i = %.10f\n", N(i), oi(i));
%    [co,~,cAb(i),errGb(i),db(i)] = stab4Poif(N(i),k,Re,'y');
%    oib(i) = imag(co)*k;
%    fprintf("N = %3d, omega_i = %.10f\n", N(i), oi(i));
end
%% Plot error compared with Orzag(1971)
oref = 0.00373967;
% oref = 0.00031866;
errn = abs(oi-oref);
% erry = abs(oib-oref);
fig1 = figure('position',[50,50,1280,720]);
h1(1) = plot(N,errn,'-ko','Displayname','Original');
hold on;
% h1(2) = plot(N,erry,'--ko','HandleVisibility','off');
% scatter(N(db==1),erry(db==1),25,'ko','filled','Displayname','Balanced');
% hold off;
set(h1,'markersize',5,'linewidth',1);
set(gca, 'YScale', 'log');
% % title('$\omega_{i,ref} = 0.00373967,\quad Re = 1\times 10^4,\quad k = 1$','fontsize',16,'interpreter','latex');
% title('$\omega_{i,ref} = 0.00031866,\quad Re = 1\times 10^9,\quad k = 0.2$','fontsize',16,'interpreter','latex');
xlabel('N','fontsize',14);
xticks([20 100:100:N(end)]);
xlim([min(N) max(N)]);
ylabel('$\left| \omega_i-\omega_{i,ref} \right|$','interpreter','latex','fontsize',14);
% legend('location','northeast','fontsize',14);
grid on;