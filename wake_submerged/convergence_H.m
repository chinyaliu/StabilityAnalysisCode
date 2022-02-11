clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,~,H,k,Fr2,Re,eps,~,~,f] = pars_wake(2);
N = 300:100:1500;
nh = linspace(0.5,5,20);
eig_spectrum = 'max';

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
        c0 = min(0.99, 1./sqrt(k(i)*Fr2));
        addvar = struct('zL1',case1.criticalH(c0),'eps',eps);
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
col = {[0 0 0],[0 0 1],[1 0 0],'#4ab712'};
lsy = {'-o','--o',':o','-.o'};
dall = abs(oih-oih(:,end));
figure;
hh = semilogy(nh, dall(1,:), lsy{1}, 'color',col{1},'Displayname',sprintf('$k = %.1f$',k(1)),'linewidth',2);
hold on;
for i = 2:lnk
    nam = sprintf('$k = %.1f$',k(i));
    semilogy(nh, dall(i,:), lsy{i}, 'Displayname', nam, 'color',col{i},'linewidth',2);
end
hold off;
xlabel('$h/\lambda$');
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$');
legend('location','northeast');
xlim([0.5 5]);
ylim([1e-16 1e-2]);
grid on;

function [oi, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi = imag(o);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
            if i~=1
                if (abs((oi-otemp)/oi)<1e-8 || i == length(N))
                    Nc = N(i);
                    break;
                end
            else
                otemp = imag(o);
            end
        end
    end
end