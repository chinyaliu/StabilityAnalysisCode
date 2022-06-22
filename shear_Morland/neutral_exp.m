clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end

%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,ud_nd,~,lamtemp,~,htemp,ff,epss,Re] = pars_Morland(1,'inv');
eig_spectrum = 'max';
hh = @(lam) lam*htemp/lamtemp;
np = 100;
delt = linspace(0,0.8,np);

%% Find analytical neutral boundaries
lmin = nan(1,np);
lmax = nan(1,np);
lp = 3;
nmin = nan(1,np);
nmax = nan(1,np);
parfor i = 1:np
    [llmin,llmax] = findneutral(ud_nd,delt(i));
    if ~isnan(llmax)
        if (strcmpi(bflow,'error function'))
            llmax = llmax+0.3;
        end
        lmin(i) = llmin;
        lmax(i) = llmax;
        lrange = linspace(llmin,llmax,lp+2);
        lpts = lrange(2:end-1);
        og = zeros(1,lp);
        for j = 1:lp
            flow1 = wMorland(N,hh(lpts(j)),ud_nd,delt(i),lpts(j),method,bflow,Re);
            c0 = min(sqrt(0.5*(lpts(j)+1./lpts(j))), 0.99*ud_nd);
            addvar = struct('zL1',flow1.invbf(c0),'eps',epss);
            o = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, ff, addvar);
            if imag(o)>0
                og(j) = 1;
            end
        end
        if sum(og) == 0
            lmin(i) = nan;
        else
            nmin(i) = lpts(find(og,1,'first'));
            nmax(i) = lpts(end-find(flip(og),1,'first')+1);
        end
    end
end
ind = ~isnan(lmin);
lmin = lmin(ind);
lmax = lmax(ind);
nmin = nmin(ind);
nmax = nmax(ind);
delt = delt(ind);
np2 = length(delt);

%%
% netmin = nan(1,np2);
% netmax = nan(1,np2);
% parfor i = 1:np2
%     flow1 = wMorland(N,1,ud_nd,delt(i),1,method,bflow,Re);
%     netmin(i) = bisearchnan(lmin(i),nmin(i),@solveflow,flow1,alg,de_singularize,do_balancing,eig_spectrum,ff,epss);
%     netmax(i) = bisearchnan(lmax(i),nmax(i),@solveflow,flow1,alg,de_singularize,do_balancing,eig_spectrum,ff,epss);
% end

%% Method 2
netmd = nan(2,np2);
bm(1) = lmin(1); bm(2) = lmax(1);
tm(1) = nmin(1); tm(2) = nmax(1);
dlam = [0.05 -0.05];
parfor j = 1:2
    for i = 1:np2
        flo(j) = wMorland(N,1,ud_nd,delt(i),1,method,bflow,Re);
        netmd(j,i) = bisearchnan(bm(j),tm(j),@solveflow,flo(j),alg,de_singularize,do_balancing,eig_spectrum,ff,epss);
        if i ~= np2
            flo(j).setprop('lambda',netmd(j,i),'h',2*netmd(j,i),'delta',delt(i+1));
            c0 = min(sqrt(0.5*(netmd(j,i)+1./netmd(j,i))), 0.99*ud_nd);
            addvar = struct('zL1',flo(j).invbf(c0),'eps',epss);
            o = flo(j).solver(alg, de_singularize, do_balancing, eig_spectrum, ff, addvar);    
            if isnan(o)
                bm(j) = netmd(j,i);
                newlam = bm(j);
                olam = nan;
                for pp = 1:10
%                 while (isnan(olam) && newlam>0 && newlam<3)
                    if (isnan(olam) && newlam>0 && newlam<3)
                    newlam = newlam+dlam(j);
                    flo(j).setprop('lambda',newlam,'h',2*newlam,'delta',delt(i+1));
                    c0 = min(sqrt(0.5*(newlam+1./newlam)), 0.99*ud_nd);
                    addvar = struct('zL1',flo(j).invbf(c0),'eps',epss);
                    olam = flo(j).solver(alg, de_singularize, do_balancing, eig_spectrum, ff, addvar);
                    else   
                        break;
                    end
                end
                tm(j) = newlam;
            else
                tm(j) = netmd(j,i);
                newlam = tm(j);
                olam = o;
                for pp = 1:10
%                 while (~isnan(olam) && newlam>0 && newlam<3)
                    if (~isnan(olam) && newlam>0 && newlam<3)
                    newlam = newlam-dlam(j);
                    flo(j).setprop('lambda',newlam,'h',2*newlam,'delta',delt(i+1));
                    c0 = min(sqrt(0.5*(newlam+1./newlam)), 0.99*ud_nd);
                    addvar = struct('zL1',flo(j).invbf(c0),'eps',epss);
                    olam = flo(j).solver(alg, de_singularize, do_balancing, eig_spectrum, ff, addvar);    
                    else
                        break;
                    end
                end
                bm(j) = newlam;
            end
        end
    end
end

%% Plot
delt1 = [linspace(0,0.18,500) linspace(0.18,delt(end),700)];
for i = 1:length(delt1)
    ldt(i) = findneutral(ud_nd,delt1(i));
end
fig1 = pltneutral(ud_nd,delt1);
hold on;
% figure;
% plot(delt1,ldt,'-k.');
% plot(delt,lmin,'-k.');
% plot(delt,lmax,'-k.');
for i = 1:2
    plot(delt,netmd(i,:),'-r.');
end
hold off; box on;
xlabel('$\Delta$');
ylabel('$\lambda$','rotation',0, 'HorizontalAlignment','right');

%%
function [oi,zc] = solveflow(lam,flow1,alg, de_singularize, do_balancing, eig_spectrum, f, epss,varargin)
    if isempty(varargin)
        c0 = min(sqrt(0.5*(lam+1./lam)), 0.99*flow1.ud);
        addvar = struct('zL1',flow1.invbf(c0),'eps',epss);
    else
        addvar = struct('zL1',varargin{1},'eps',epss);
    end
    flow1.setprop('lambda',lam,'h',2*lam);
    o = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    if isnan(o)
        oi = nan;
    else
        oi = imag(o(1));
    end
    zc = flow1.zc;
end