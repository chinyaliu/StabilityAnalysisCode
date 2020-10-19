close all; clear all;% clc
load('neutral_pts.mat');
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 401;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
if (strcmpi(method(1),'d4') && strcmpi(method(2),'schimd') && strcmpi(method(3),'d4'))
    M = N-2;
else
    M = N;
end
z = cospi((0:1:M)/M)'; 
U = 1-z.^2;

R_tar = [1e4 1e5 1e6 1e7 1e8 1e9];
for i = 1:length(R_tar)
    [~,pos] = min(abs(Re0-R_tar(i)));
    if ( (pos+5)>length(Re0) || (pos-5)<1)
        k_tar = sort(k0(pos-10:pos));
    else
        k_tar = sort(k0(pos-5:pos+5));
    end
    kn = k_tar(1)-0.01:0.01:k_tar(end)+0.01;
    for j = 1:length(kn)
        [o,~,~,~,~] = poiseuille_solver(N,kn(j),R_tar(i),method,alg,do_balancing);
        R(i).oi(j) = imag(o);
        R(i).or(j) = real(o);
        if (imag(o)<0)
            R(i).or(j) = nan;
        end
    end
    R(i).kn = kn;
end

fig1 = figure('position',[50,50,1280,1440]);
subplot(2,1,1);
for i = 1:length(R_tar)
    text = sprintf('%.0e',R_tar(i));
    plot(R(i).kn,R(i).or,'marker','.','markersize',10,'Displayname',text,'linewidth',1);
    hold on;
end
hold off;
legend('location','northwest','fontsize',24);
set(gca,'fontsize',24);
% xlabel('k','fontsize',30);
ylabel('\omega_r ','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
grid on;
subplot(2,1,2);
for i = 1:length(R_tar)
    text = sprintf('%.0e',R_tar(i));
    plot(R(i).kn,R(i).oi,'marker','.','markersize',10,'Displayname',text,'linewidth',1);
    hold on;
end
hold off;
% legend('location','northwest','fontsize',24);
set(gca,'fontsize',24);
ylim([0 0.008]);
xlabel('k','fontsize',30);
ylabel('\omega_i ','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -3;
grid on;

fig2 = figure('position',[50,50,1280,720]);
for i = 1:length(R_tar)
    text = sprintf('%.0e',R_tar(i));
    c_r = R(i).or./R(i).kn;
    plot(R(i).kn,c_r,'marker','.','markersize',10,'Displayname',text,'linewidth',1);
    hold on;
end
hold off;
legend('location','southeast','fontsize',24);
set(gca,'fontsize',24);
xlabel('k','fontsize',30);
ylabel('c_r ','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
grid on;