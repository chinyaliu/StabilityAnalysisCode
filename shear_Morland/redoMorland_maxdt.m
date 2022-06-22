clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    

%% Inputs
[method,alg,bflow,de_singularize,do_balancing,~,N,~,~,~,~,~,f,epss,Re] = pars_Morland(1,'inv');
eig_spectrum = 'max';

%% Morland's data ( as input reference )
if strcmp(bflow,"exponential")
    ud_nd = [1.75 2 2.25 2.5];
    delta_nd = [0.372 0.291 0.235 0.194];
    lam_nd = [0.870 0.817 0.767 0.722]; 
elseif strcmp(bflow,"error function")
    ud_nd = [1.75 2 2.25 2.5];
    delta_nd = [0.299 0.234 0.187 0.150];
    lam_nd = [0.823 0.764 0.707 0.649]; 
end

%% Find maximum for each U_0
tic;
ddt = 0.01;
dlam = 0.01;
sol_var = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',-0.1,'eps',epss)};
dt_list = NaN(1,length(ud_nd));
lam_list = NaN(1,length(ud_nd));
o_list = NaN(1,length(ud_nd));
parfor i = 1:length(ud_nd)
    flow1 = wMorland(N,2*lam_nd(i),ud_nd(i),delta_nd(i),lam_nd(i),method,bflow,Re);
    oref = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',-0.1,'eps',epss));
    advar = {lam_nd(i),dlam,flow1,alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',flow1.zc,'eps',epss)};
    delta = goldensearchmax(delta_nd(i)-ddt,delta_nd(i)+ddt,@maxlambda,advar{:});
    if ~isnan(delta)
        flow1.setprop('delta',delta);
        lambda = goldensearchmax(lam_nd(i)-dlam,lam_nd(i)+dlam,@solveflow,advar{3:end});
        if ~isnan(lambda)
            flow1.setprop('lambda',lambda,'h',2*lambda);
            o = flow1.solver(advar{4:end});
            o_list(i) = o(1);
            dt_list(i) = delta;
            lam_list(i) = lambda;
        end
    end
end
c_list = 0.5*o_list.*lam_list/pi;
toc;

function oimax = maxlambda(delta,lam,dlam,flow1,alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    flow1.setprop('delta',delta);
    xL = lam-dlam;
    xU = lam+dlam;
    if ~isnan(xL)
        advar = {flow1,alg, de_singularize, do_balancing, eig_spectrum, f, addvar};
        lambda = goldensearchmax(xL,xU,@solveflow,advar{:});
        flow1.setprop('lambda',lambda,'h',2*lambda);
        oo = flow1.solver(advar{2:end});
        oimax = imag(oo(1));
    end
end

function oi = solveflow(lam,flow1,alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    flow1.setprop('lambda',lam,'h',2*lam);
    o = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    oi = imag(o(1));
end