clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,H,k,Fr2,Re,eps,c0,~,f] = pars_wake(2);
N = 300:100:1500;
nh = linspace(0.5,5,20);
%% Run solver
lnk = length(k);
lnh = length(nh);
oih = nan(lnk,lnh);
Nc = nan(lnk,lnh);
tic;
parfor i = 1:lnk
    for j = 1:lnh
        h = nh(j)*2*pi./k(i)+H;
        case1 = wSubmerged(100,H,k(i),h,Re,Fr2);
        case1.numMeth(method);
        addvar = struct('zL1',0.74708299+H,'eps',eps);
        [oih(i,j), Nc(i,j)] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    end
end
toc;
% %% Plot diff(c_i) vs h
% dc = abs(diff(imag(oih)));
% figure;
% semilogy(nh(1:end-1),dc,'-bo');
% xlabel('$n\lambda$');
% ylabel('$\| \omega_i(m)-\omega_i(m+1) \|$');
% grid on;
%% Plot figure
col = {[0 0 0],[0 0 1],[1 0 0],"#77AC30"};
dall = abs(oih-oih(:,end));
figure;
h = semilogy(nh, dall(1,:), '-ko','Displayname',sprintf('$k = %.2f$',k(1)));
hold on;
for i = 2:lnk
    nam = sprintf('$k = %.2f$',k(i));
    semilogy(nh, dall(i,:), '-o', 'Displayname', nam, 'color',col{i});
end
hold off;
xlabel('$h/\lambda$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |$');
legend('location','northeast');
xlim([0.5 5]);
ylim([1e-16 1e-2]);
grid on;
% %% Plot figure
% figure;
% nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oih(end))];
% semilogy(nh, oih, '-o');
% xlabel('$h/\lambda$');
% ylabel('$\omega_i$');
% nam = [sprintf('$k = %.2f, ',k) '\' sprintf('omega = %0.9e$',oih(end))];
% legend(nam,'location','northeast');
% grid on;
% %% Plot figure
% figure;
% nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oih(end))];
% semilogy(nh, Nc, '-o');
% xlabel('$h/\lambda$');
% ylabel('grids required to converge');
% nam = [sprintf('$k = %.2f, ',k) '\' sprintf('omega = %0.9e$',oih(end))];
% legend(nam,'location','southeast');
% grid on;

function [oi, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi = imag(o);
        if i~=1
            if (abs((oi-otemp)/oi)<1e-10 || i == length(N))
                Nc = N(i);
                break;
            end
        else
            otemp = imag(o);
        end
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
        end
    end
end