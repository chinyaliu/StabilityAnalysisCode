clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end

%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,~,~,~,~,~,ff,epss,~] = pars_Morland(1,'vis');
eig_spectrum = 'all';
addvar = struct('zL1',-0.1,'eps',epss);
sol_var = {alg, de_singularize, do_balancing, eig_spectrum, ff, addvar};

%% Set search domain and method
lamlim = [0 20];
ulim = [0 100];
m = [5 10];
tol = 1e-5;

%% 
lam = cell(2,1);
u = cell(2,1);
parfor j = 1:2
    md = choosemode(j);
    % find three base points
    mb = [0 m(j) inf];
    base0 = [ulim(end) 2 lamlim(end)];
    u_ini = nan(1,length(mb));
    lam_ini = nan(1,length(mb));
    for i = 1:length(mb)
        p = setline(mb(i),base0(i));
        sol = sol_var;
        mdmode = md;
        pbd = [1 20];
        if (i~=1) && (j==2)
            pbd = [10 20];
        end
        flow1 = wMorland(N,1,1,1,1,method,bflow,1);
        p_ini = secantf(pbd(1),pbd(2),p,flow1,sol,mdmode,tol);
        [u_ini(i),lam_ini(i)] = p(p_ini);
    end
    % four workers
    p_ini = invp(m(j),u_ini,lam_ini);
    ntotal = 100;
    dlam = 2*diff(p_ini)/ntotal;
    dp = [dlam(1) dlam(1) dlam(2) dlam(2)];
    dlam0 = [dlam(1) -dlam(1) -dlam(2) dlam(2)];
    lam0 = [p_ini p_ini(2)];
    p0 = [lam_ini lam_ini(2)];
    u_list = nan(4,ntotal/4);
    lam_list = nan(4,ntotal/4);
    for i = 1:4
        flow1 = wMorland(N,1,1,1,1,method,bflow,1);
        [ui,lami] = solset(flow1,m(j),lam0(i),dlam0(i),p0(i),dp(i),sol_var,md,tol,ntotal);
        u_list(i,:) = ui;
        lam_list(i,:) = lami;
    end
    utemp = [u_list(:); u_ini.'];
    lamtemp = [lam_list(:); lam_ini.'];
    [lam{j},ind] = sort(lamtemp);
    u{j} = utemp(ind);
end

%%
figure('position',[50,0,1080,810]);
plot(lam{1},u{1},'b');
hold on; plot(lam{2},u{2},'r'); hold off;
% area(lam{1},u{1},ulim(end),'FaceColor','b');
% hold on; area(lam{2},u{2},ulim(end),'FaceColor','r'); hold off;
grid on;
xlim([0 20]); 
ylim([0 100]);
axis square;
xlabel('$\lambda$(cm)');
ylabel('$u_{*}$(cm/s)');
legend('mode 1','mode 2','location','east');

%%
function [u,lam] = solset(flow1,m,lam0,dlam0,p_in,dp,sol,mdmode,tol,ntotal)
    u = nan(1,ntotal/4);
    lam = nan(1,ntotal/4);
    p = p_in;
    for i = 1:ntotal/4
        pxy = setline(m,lam0+i*dlam0);
        pp = [p-dp p+dp];
        pp(pp<0) = p;
        p_out = secantf(pp(1),pp(2),pxy,flow1,sol,mdmode,tol);
        if isnan(p_out)
            pp = [1 20];
            p_out = secantf(pp(1),pp(2),pxy,flow1,sol,mdmode,tol);
        end
        [u(i),lam(i)] = pxy(p_out);
        p = 2*p_out-p;
    end
end

% secant method
function p_out = secantf(p0,p1,pxy,flow1,sol_var,md,tol)
    % initialize
    p_out = nan;
    f0 = findf(p0);
    f1 = findf(p1);
    if f0<f1
        % iterations
        for i = 1:500 % maximum iteration
            p2 = p1 - f1*(p0-p1)/(f0-f1);
            f2 = findf(p2);
            if abs(f2) < tol
                p_out = p2;
                break;
            end
%             if (f2>0) && (f2<f1)
%                 p1 = p2;
%                 f1 = f2;
%             elseif (f2<0) && (f2>f0)
%                 p0 = p2;
%                 f0 = f2;
%             else
%                 error('New eigenvalue out of bounds');
%             end
            if (f2<0) && (f2>f0)
                p0 = p2;
                f0 = f2;
            else
                p1 = p2;
                f1 = f2;
            end
        end
    end

    function f = findf(p)
        [u,lam] = pxy(p);
        [ud,del,lamd,Re] = wtoMor(u,lam);
        h = lamd;
        flow1.setprop('ud',ud,'delta',del,'lambda',lamd,'Re',Re,'h',h);
        o = flow1.solver(sol_var{:});
        f = md(o./flow1.k,flow1.ud);
    end
end

% search line
function setxy = setline(m,lam0)
%
% Dimensions u*(cm/s) lambda(cm)
% m   : slope 
% lam : lambda when u* = 0
%
    if m == 0
        setxy = @p0;
    elseif isinf(m)
        setxy = @pinf;
    else
        setxy = @p;
    end
    function [u, lam] = p(lam_in)
        u = m*(lam_in-lam0);
        lam = lam_in;
    end
    function [u,lam] = p0(lam_in)
        u = lam0;
        lam = lam_in;
    end
    function [u,lam] = pinf(lam_in)
        u = 5*lam_in;
        lam = lam0;
    end
end

function lam0 = invp(m,u,lam)
    lam0 = lam-u./m;
end

% choose mode
function md = choosemode(chm)
    switch(chm)
        case(1)
            md = @md1;
        case(2)
            md = @md2;
    end
    function ci = md1(c,varargin)
        c_neg = c(real(c)<0);
        if ~isempty(c_neg)
            [~, ind] = max(imag(c_neg));
            ci = imag(c_neg(ind));
        else
            error('Cannot determine mode 1');
        end
    end
    function ci = md2(c,ud_nd)
        % find discrete eigenvalues
        a = 1:length(c); 
        crange = ((real(c)>-1e-5) & (real(c)-ud_nd<=1e-5)) & imag(c) > -1;
        aa = a(crange);
        abci = isoutlier(imag(c(aa)),'movmedian',5);
        abcr = isoutlier(real(c(aa)),'movmedian',5);
        aa = aa(abci|abcr);
        c_chosen = c(aa);
        [~, ind] = max(real(c_chosen));
        ci = imag(c_chosen(ind));
    end
end