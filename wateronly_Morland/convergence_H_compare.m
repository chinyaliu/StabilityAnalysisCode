clear all; close all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,~,bflow,de_singularize,do_balancing,eig_spectrum,~,ud_nd,delta_nd,lambda_nd,c0,~,f] = pars_Morland(2);
alg = "qr";
alg2 = "invB";
alg3 = "qz";
% init_var = {1000,0,ud_nd,delta_nd,lambda_nd,method,bflow};
% sol_var = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};
% sol_var2 = {alg2, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};
% sol_var3 = {alg3, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};

%% Run solver
h_list = linspace(0.5,5,20);
N_list = 300:100:1500;
% N_list = 800;
c_list = NaN(1,length(h_list));
c_list2 = NaN(1,length(h_list));
c_list3 = NaN(1,length(h_list));
h_Nlist = NaN(1,length(h_list));
h_Nlist2 = NaN(1,length(h_list));
h_Nlist3 = NaN(1,length(h_list));
tic;
parfor i = 1:length(h_list)
    case1 = wMorland(1000,0,ud_nd,delta_nd,lambda_nd,method,bflow);
    case1.h = h_list(i)*lambda_nd;
    
    [ci, hN] = convgmode(N_list, case1, alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2));
    c_list(i) = ci;
    h_Nlist(i) = hN;
    
    [ci, hN] = convgmode(N_list, case1, alg2, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2));
    c_list2(i) = ci;
    h_Nlist2(i) = hN;

    [ci, hN] = convgmode(N_list, case1, alg3, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2));
    c_list3(i) = ci;
    h_Nlist3(i) = hN;
end
toc;

%% Plot N vs h
figure;
plot(h_list,h_Nlist,'-bo');
hold on;
plot(h_list,h_Nlist2,'-ro');
plot(h_list,h_Nlist3,'-go');
hold off;
xlabel('$n\lambda$');
ylabel('$N$');
legend('eig(A\textbackslash B)','eig(B\textbackslash A)','qz(A,B)','location','northwest');
grid on;

%% Plot c_i vs h
figure;
plot(h_list,c_list,'-bo');
hold on;
plot(h_list,c_list2,'-ro');
plot(h_list,c_list3,'-go');
hold off;
xlabel('$n\lambda$');
ylabel('$c_i$');
legend('eig(A\textbackslash B)','eig(B\textbackslash A)','qz(A,B)','location','northwest');
grid on;

%% Plot diff(c_i) vs h
dc = abs(diff(imag(c_list)));
dc2 = abs(diff(imag(c_list2)));
dc3 = abs(diff(imag(c_list3)));
figure;
semilogy(h_list(1:end-1),dc,'-bo');
hold on;
semilogy(h_list(1:end-1),dc2,'-ro');
semilogy(h_list(1:end-1),dc3,'-go');
hold off;
xlabel('$n\lambda$');
ylabel('$\| c_i(m)-c_i(m+1) \|$');
legend('eig(A\textbackslash B)','eig(B\textbackslash A)','qz(A,B)','location','northwest');
grid on;

function [ci, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    ctemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solvers(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        ci = imag(o);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
            if i~=1
                if (abs((ci-ctemp)/ci)<1e-8 || i == length(N))
                    Nc = N(i);
                    break;
                else
                    ctemp = ci;
                end
            else
                ctemp = imag(o);
            end
        end
    end
end