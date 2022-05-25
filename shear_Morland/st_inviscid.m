clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end

%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,N,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(3);
eig_spectrum = 'all';
% Cusp method from k_i
ki = [-0.1 -0.5 -1 -1.5];

%% Real k
tic;
p1 = wMorland(N,h(1),ud_nd,delta_nd,lambda_nd(1),method,bflow,Re);
addvar = struct('zL1',p1.invbf(c0(1)),'eps',epss);
or = []; kr = [];
cutz = p1.invbf(c0);
cut_temp = cutz(1);
for i = 1:length(lambda_nd)
    p1.setprop('lambda',lambda_nd(i),'h',h(i));
    oall = p1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    or = [or; oall]; 
    kr = [kr; (2*pi./lambda_nd(i))*ones(length(call),1)];
    if isnan(p1.zc)
        cutz(i) = cut_temp;
        addvar.zL1 = cut_temp;
    else
        cut_temp = p1.zc;
        cutz(i) = cut_temp;
        addvar.zL1 = cut_temp;
    end
end
cr = or./kr;

% %% Run solver
% tic;
% cc = cell(length(ki),1);
% kk = cell(length(ki),1);
% cutt = cutz(1);
% parfor j = 1:length(ki)
%     addvar = struct('zL1',cutt,'eps',eps);
%     p1 = wMorland(N,2,ud_nd,delta_nd,2,method,bflow);
%     [c1,o1] = findmodes(2*pi./lambda_nd+ki(j)*1i,h,addvar,p1,alg, de_singularize, do_balancing, eig_spectrum, f,cutz);
%     cc{j} = c1;
%     kk{j} = o1;
% end
% toc;

%% Run solver 2
tic;
cc = cell(length(ki),1);
kk = cell(length(ki),1);
parfor j = 1:length(ki)
    cutt = cutz; ht = h;
    for i = 1:length(lambda_nd)
        addvar = struct('zL1',cutt(i),'eps',epss);
        p1 = wMorland(N,ht(i),ud_nd,delta_nd,2,method,bflow,Re);
        p1.setprop('k',2*pi./lambda_nd(i)+ki(j)*1i);
        c = p1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        cc{j} = [cc{j}; c];
        kk{j} = [kk{j}; (2*pi./lambda_nd(i)+ki(j)*1i)*ones(length(c),1)];
    end
end
toc;

%%
figure; 
ax = axes();
plot(real(cr),imag(cr),'k.','Markersize',8,'HandleVisibility','off');
hold on;
for i = 1:length(ki)
    hh(i) = plot(real(cc{i}),imag(cc{i}),'.','Markersize',8,'DisplayName',sprintf('$%+.2f$',ki(i)));
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
% ylim([-0.5 0.3]);
% xlim([-1 2]);
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% Set legend
hCopy = copyobj(hh, ax); 
for i = 1:length(ki)
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 20; 
end
leg = legend(hCopy,'location','southeastoutside');
title(leg,'$k_i$');

%% o
figure; 
ax = axes();
plot(real(or),imag(or),'k.','Markersize',8,'HandleVisibility','off');
hold on;
for i = 1:length(ki)
    oo{i} = cc{i}.*kk{i};
    hh(i) = plot(real(oo{i}),imag(oo{i}),'.','Markersize',8,'DisplayName',sprintf('$%+.2f$',ki(i)));
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
% ylim([-0.5 0.3]);
% xlim([-1 2]);
xlabel('$\omega_r$','fontsize',30);
ylabel('$\omega_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% Set legend
hCopy = copyobj(hh, ax); 
for i = 1:length(ki)
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 20; 
end
leg = legend(hCopy,'location','southeastoutside');
title(leg,'$k_i$');
    
function [cc,kc] = findmodes(kk,hh,addvar,p1,alg, de_singularize, do_balancing, eig_spectrum, f)
    cc = []; kc = [];
    for i = 1:length(kk)
        p1.setprop('k',kk(i),'h',hh(i));
        c = p1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        if real(c(1)) > 0 && real(c(1)) <= p1.ud
            addvar.zL1 = p1.invbf(real(c(1)));
        end
        cc = [cc; c];
        kc = [kc; kk(i)*ones(length(c),1)];
    end
end