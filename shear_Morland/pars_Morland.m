function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lam_nd,c0,h,f,eps] = pars_Morland(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(4);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'all';
N = 600;
ud_nd = 2;
delta_nd = 0.291;
lam_nd = 0.817;
% delta_nd = 0.234;
% lam_nd = 0.764;
ddm_number = 44;
eps = 0.1;
f = wMorland.ddmtype(ddm_number);
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
%             % exponential
%             ud_nd = [1.75 2 2.25 2.5];
%             delta_nd = [0.372 0.291 0.235 0.194];
%             lam_nd = [0.870 0.817 0.767 0.722]; 
            % error function
            ud_nd = [1.75 2 2.25 2.5];
            delta_nd = [0.299 0.234 0.187 0.150];
            lam_nd = [0.823 0.764 0.707 0.649]; 
        case(3)
            lam_nd = linspace(0.01,2,300);
    end
end
h = 2*lam_nd;
c0 = min(sqrt(0.5*lam_nd), 0.99*ud_nd);

end