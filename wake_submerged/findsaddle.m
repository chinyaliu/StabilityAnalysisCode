clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[solfunc,Re] = setfunc;
fsimp = setsimplex(solfunc,Re);
% test
k = 1.6-1i;
oin = solfunc(k);
fig = figure;
plot(real(k),imag(k),'ro');
xlim([0 4]);
ylim([-2 0]);
[kmin,fval] = fminsearch(fsimp,[real(k),-imag(k)],optimset('Display','iter','MaxIter',100,'OutputFcn',@outfun));
plot(real(kmin),imag(kmin),'bo');
hold off;
oout = solfunc(kmin(1)-1i*kmin(2));
fmin = popk(solfunc,kmin);
%% Broyden's method
% k = 0.3-1.5i;
[kmin2, fmin2] = ntmeth(solfunc,[real(k),-imag(k)]);
oout2 = solfunc(kmin2);

%%
function [fout,Re] = setfunc
    [method,alg,bflow,de_singularize,do_balancing,~,N,H,k1,Fr2,Re,epss,~,h,f] = pars_wake(1,'inv');
    f = wSubmerged.ddmtype(1);
    eig_spectrum = 'all';
    flow1 = wSubmerged(N,H,k1,h,Re,Fr2,bflow);
    flow1.numMeth(method);
    if isinf(Re)
        fout = @solveRE;
    end

    function out = solveRE(k)
        c0 = min(0.99, 1./sqrt(real(k)*flow1.Fr2));
        addvar = struct('zL1',flow1.invbf(c0),'eps',epss);
        flow1.setprop('k',k,'h',2*2*pi./real(k));
        o = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        op = o(real(o)>0);
        c = op./k;
        [mc,ind] = max(imag(c));
        if mc<1e-5
%             error('Wrong eigenvalue\n');
            out = nan+nan*1i;
        else
            out = op(ind);
        end
    end
end

function fout = setsimplex(solfunc,Re)
if isinf(Re)
    fout = @fsimplex;
end
    function out = fsimplex(kin)
        k = kin(1)-1i*kin(2);
        delt = 1e-6;
        o = solfunc(k);
        or = solfunc(k+delt);
        oi = solfunc(k-delt*2i);
        out = ((imag(or)-imag(o))/delt)^2+(0.5*(imag(o)-imag(oi))/delt)^2;
        if isnan(out)
            out = 100;
        end
    end
end

function stop = outfun(x, ~, ~)
    stop = false;
    hold on;
    plot(x(1),-x(2),'kx');
    drawnow
end

function [out,fout] = popk(solfunc,kin)
    k = kin(1)-1i*kin(2);
    dtr = 1e-6;
    dti = dtr*10;
    o = solfunc(k);
    or = solfunc(k+dtr);
    oi = solfunc(k-dti*1i);
    out = 0.5*((or-o)/dtr-1i*(o-oi)/dti);
end

function out = ntmeth(solfunc,kin)
    % Broyden's method
    k0 = kin(1)-1i*kin(2);
    f0 = f(k0);
    fp0 = fp(k0);
    kn = k0-f0/fp0;
    fn = f(kn);
    fpn = (fn-f0)/(kn-k0);
    for i = 1:100
        if abs(fn)<1e-8
            break;
        end
        kn1 = kn-fn/fpn;
        fn1 = f(kn1);
        fpn = (fn1-fn)/(kn1-kn);
        fn = fn1;
        kn = kn1;
        fprintf('iteration %2d, f = %.3e, k = %.3f-%.3fi\n',i,abs(fn), real(kn), -imag(kn))
    end
    out = kn;
    fout = fn;

    function out = f(k)
        dtr = 1e-6;
        dti = dtr*10;
        o = solfunc(k);
        or = solfunc(k+dtr);
        oi = solfunc(k-dti*1i);
        out = 0.5*((or-o)/dtr-1i*(o-oi)/dti);
    end

    function out = fp(k)
        dtr = 1e-6;
        dti = dtr*10;
        o = solfunc(k);
        or = solfunc(k+dtr);
        ro = solfunc(k-dtr);
        frr = (or-2*o+ro)/dtr^2;
        oi = solfunc(k+dti*1i);
        io = solfunc(k-dti*1i);
        fii = (oi-2*o+io)/dti^2;
        ori = solfunc(k+dtr+dti*1i);
        roi = solfunc(k-dtr+dti*1i);
        ior = solfunc(k+dtr-dti*1i);
        iro = solfunc(k-dtr-dti*1i);
        fri = (ori-roi-ior+iro)/4/dtr/dti;
        out = 0.25*(frr-fii-2i*fri);
    end
end