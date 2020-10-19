close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 400;
% k = 1;
% Re = 1e4;
k = 0.2;
Re = 1e9;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
[o, an, cA, errGEP, db] = poiseuille_solver(N,k,Re,method,alg,do_balancing);
%% Find critical height
if (strcmpi(method(1),'d4') && strcmpi(method(2),'schimd') && strcmpi(method(3),'d4'))
    M = N-2;
else
    M = N;
end
[z,D] = Dcheb(N,1,method(1));
U = 1-z.^2;
[minpt,minpos] = min(abs(U(1:fix(N+1)/2)-real(o)/k));
%% Get streamfunction phi
switch lower(method(1))
    case {'d2', 'd4'}
        phi = D(:,:,1)*an(1:N+1);
        up = D(:,:,2)*an(1:N+1);
        wp = -1i*k*phi;
    case 'uw'
        up = D0*an(1:N+1);
        wp = D0*an(N+2:end);
        phi = -wp/k/(1i);
end
phi = phi./max(abs(phi));
up = up./max(abs(phi));
wp = wp./max(abs(phi));
%% Plot phi
fig1 = figure('position',[0 0 1260 640]);
subplot(1,4,1);
plot(real(phi(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(real(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('real(\phi)','FontSize',22);
ylabel('z','FontSize',22);
grid on;
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
subplot(1,4,2);
plot(imag(phi(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(imag(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(\phi)','FontSize',22);
ylabel('z','FontSize',22);
% title('imag','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,3);
plot(abs(phi(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(abs(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('abs(\phi)','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,4);
plot(angle(phi(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(angle(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('angle(\phi)','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
sgtitle('\phi','FontSize',22);
saveas(fig1,'phi.png');
%% Plot u
fig2 = figure('position',[0 0 1260 640]);
subplot(1,4,1);
plot(real(up(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(real(up(minpos)),z(minpos),'ro'); hold off;
xlabel('real(u'')','FontSize',22);
ylabel('z','FontSize',22);
% title('real','FontSize',22);
grid on;
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
subplot(1,4,2);
plot(imag(up(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(imag(up(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(u'')','FontSize',22);
ylabel('z','FontSize',22);
% title('imag','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,3);
plot(abs(up(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(abs(up(minpos)),z(minpos),'ro'); hold off;
xlabel('abs(u'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,4);
plot(angle(up(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(angle(up(minpos)),z(minpos),'ro'); hold off;
xlabel('angle(u'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
sgtitle('u''','FontSize',22);
saveas(fig2,'u_per.png');
%% Plot w
fig3 = figure('position',[0 0 1260 640]);
subplot(1,4,1);
plot(real(wp(1:fix(N+1)/2)),z(1:fix(N+1)/2),'k','linewidth',1);
hold on; plot(real(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('real(w'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,2);
plot(imag(wp(1:fix(N+1)/2)),z(1:fix(N+1)/2),'k','linewidth',1);
hold on; plot(imag(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(w'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,3);
plot(abs(wp(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(abs(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('abs(w'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
subplot(1,4,4);
plot(angle(wp(1:fix(N+1)/2)),z(1:fix(N+1)/2),'-k','linewidth',1);
hold on; plot(angle(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('angle(w'')','FontSize',22);
ylabel('z','FontSize',22);
set(gca,'fontsize',16);
if (Re > 1e7) ylim([0.9 1]); yticks(0.9:0.02:1); end
grid on;
sgtitle('w''','FontSize',22);
saveas(fig3,'w_per.png');