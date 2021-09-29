clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,~,~,h,f] = pars_Morland(1);

%% Run solver
col = repmat(get(gca,'colororder'),5,1);
[lmin, lmax] = findneutral(ud_nd,delta_nd); % for exponential profile
lambda_list = linspace(lmax,lmin,200);
ki_list = [0 0.1 0.5 1 1.5];
o_list = cell(length(ki_list),1);
t1 = tic;

flow1 = wMorland(N,h,ud_nd,delta_nd,1,method,bflow);
addvar = struct('zL1',delta_nd,'eps',0.2);
for i = 1:length(ki_list)
    ki = ki_list(i);
    addvar.zL1 = delta_nd;
    % Select mode to observe
    cs = [4.1457;-0.35307];
    o_list{i} = [];
    for j = 1:length(lambda_list)
        lam = lambda_list(j);
        fprintf('wavelength = %.2f\n', lam);
        flow1.k = 2*pi/lam-ki*1i;
        flow1.h = 3*lam;
%         c_all = flow1.solver(alg, do_balancing, eig_spectrum, f, addvar);
        c_all = flow1.solvers(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        if isnan(flow1.zc)
            addvar.zL1 = delta_nd;
        else
            addvar.zL1 = -flow1.zc;
        end

        % Choose selected eigenvalues
        c_new = finddiscrete(c_all,ud_nd);
        c0 = repmat(cs.',length(c_new),1);
        c1 = repmat(c_new,1,length(cs));
        [~,min_ind] = min(abs(c0-c1),[],1);
        [~,min2] = min(abs(c0-c1),[],2);
        for m = 1:length(cs)
            if min2(min_ind(m)) == m
                cs(m) = c_new(min_ind(m));
            else
                cs(m) = NaN;
            end
        end
        cs(isnan(cs)) = [];
        if length(cs)<length(c_new)
            cs = [cs;c_new(setdiff(1:end,min_ind))];
        end
        o_list{i} = [o_list{i};cs*flow1.k];
    end
end
hold off;
toc(t1);

%%
figure;
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold on;
xline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
for i = 1:length(o_list)
    o = o_list{i};
    textk = sprintf('$k_i=%+1.1f$',-ki_list(i));
    plot(real(o),imag(o),'.','DisplayName',textk);
end
hold off;
xlabel('$\omega_r$');
ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
[~, objh] = legend('interpreter','latex','location','southeastoutside');
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 20);
% ylim([-0.5 0.5]);
% xlim([-5 20]);

function c = finddiscrete(ev,ulim)
    a = 1:length(ev); 
    crange = (real(ev)>0) & (real(ev)<ulim);
    aa = a(crange);
    abch = isoutlier(imag(ev(aa)),'movmedian',5);
    newa = aa(abch);
    newa(abs(imag(ev(newa)))<1e-4) = [];
    aa = [a(real(ev)<0) newa a(real(ev)>ulim)];
    c = ev(aa);
end