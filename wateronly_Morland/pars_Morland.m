function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lam_nd,c0,h,f] = pars_Morland(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(4);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(2);
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'all';
N = 600;
ud_nd = 2;
delta_nd = 0.291;
lam_nd = 0.817;
h =3*lam_nd;
ddm_number = 92;
c0 = sqrt(0.5*(lam_nd+1./lam_nd));
f = wMorland.ddmtype(ddm_number);
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
            eig_spectrum = 'max';
    end
end
end