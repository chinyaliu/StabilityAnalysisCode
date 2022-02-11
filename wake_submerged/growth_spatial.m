clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 

%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,N,H,k,Fr2,Re,eps,~,h,f] = pars_wake(3);
eig_spectrum = 'all';
% Cusp method from k_i
ki = [-0.5 -1 -1.5 -1.7 -2 -2.5 -2.7];
inflec_pt = -0.74708299-H;

%% Real k                            
tic;
addvar = struct('zL1',-inflec_pt,'eps',eps);
p1 = wSubmerged(N,H,k(1),h(1),Re,Fr2);
p1.numMeth(method);
% or = findmodes(k,h,addvar,p1,alg, de_singularize, do_balancing, eig_spectrum, f);
or = ones(length(k),1);
cutz = -inflec_pt;
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    oall = p1.solver(alg, de_singularize, do_balancing, 'max', f, addvar);
    or(i) = oall(1);
    if isnan(p1.zc)
        addvar.zL1 = cutz;
    else
        cutz = p1.zc;
        addvar.zL1 = p1.zc;
    end
end

%% Run solver
tic;
oo = cell(length(ki),1);
parfor j = 1:length(ki)
%     kk = k+ki(j)*1i;
    addvar = struct('zL1',-inflec_pt,'eps',eps);
%     hh = h;
    p1 = wSubmerged(N,H,0.01,6,Re,Fr2);
    p1.numMeth(method);
    o3 = findmodes(k+ki(j)*1i,h,addvar,p1,alg, de_singularize, do_balancing, eig_spectrum, f);
    oo{j} = o3;
end
toc;

%%
figure; 
ax = axes();
plot(real(or),imag(or),'k.','Markersize',8,'HandleVisibility','off');
hold on;
for i = 1:length(ki)
    hh(i) = plot(real(oo{i}),imag(oo{i}),'.','Markersize',8,'DisplayName',sprintf('$%+.1f$',ki(i)));
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
ylim([-0.5 0.3]);
xlim([-1 2]);
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% Set legend
hCopy = copyobj(hh, ax); 
for i = 1:length(ki)
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 20; 
end
leg = legend(hCopy,'location','southeastoutside');
title(leg,'$k_i$');
    
function o3 = findmodes(kk,hh,addvar,p1,alg, de_singularize, do_balancing, eig_spectrum, f)
    o3 = [];
    for i = 1:length(kk)
        p1.k = kk(i); p1.h = hh(i);
        oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oall = oall(real(oall)>-50); % remove artificial eigenvalues
        o = oall(1);
        if real(o) > 0
            addvar.zL1=p1.criticalH(real(o)/kk(i));
        end
        
        ca = oall/kk(i);
        % Criteria 2
        a = 1:length(ca); 
        crange = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-5));
        dis = ((real(ca)-1).^2 +imag(ca).^2)>1e-5;
        aa = a(crange&dis);
        abch = isoutlier(imag(ca(aa)),'movmedian',20);
        aa1 = [a(dis&~crange) aa(abch)];
        
        %
        creal = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-5));
        cimag = abs(imag(ca)) > 0.005;
        aa2 = a(cimag|~creal);
        
        %
        A = [aa1 aa2];
        [~,IA,~] = unique(A);
        aa3 = unique(A(setdiff((1:length(A)),IA)));
        o3 = [o3; oall(aa3)];
    end
end