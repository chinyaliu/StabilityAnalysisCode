clear all; close all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,~,bflow,de_singularize,do_balancing,eig_spectrum,~,ud_nd,delta_nd,lambda_nd,c0,~,f] = pars_Morland(2);
alg = "qr";
alg2 = "invB";
alg3 = "qz";
init_var = {1000,0,ud_nd,delta_nd,lambda_nd,method,bflow};
sol_var = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};
sol_var2 = {alg2, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};
sol_var3 = {alg3, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};

%% Run solver
h_list = linspace(0.5,5,20);
% N_list = 300:100:1500;
N_list = 800;
c_list = NaN(1,length(h_list));
c_list2 = NaN(1,length(h_list));
c_list3 = NaN(1,length(h_list));
h_Nlist = NaN(1,length(h_list));
h_Nlist2 = NaN(1,length(h_list));
h_Nlist3 = NaN(1,length(h_list));
cA1 = NaN(1,length(h_list));
cA2 = NaN(1,length(h_list));
cA3 = NaN(1,length(h_list));
tic;
case1 = wMorland(init_var{:});
for i = 1:length(h_list)
    fprintf('h = %.2f times wave length.\n',h_list(i));
    case1.h = h_list(i)*lambda_nd;
%     c_temp = 0;
%     for N = N_list
%         fprintf('N = %3d\n', N);
%         case1.N = N;
%         c = case1.solvers(sol_var{:});
%         if abs(c-c_temp)<1e-8
%             break;
%         elseif ~isnan(c_temp)
%             c_temp = c;
%         end
%         if ~isnan(case1.zc)
%             sol_var{6}.zL1 = -case1.zc;
%         end
%     end
    case1.N = N_list(end);
    [c,~,cA] = case1.solvers(sol_var{:});
    cA1(i) = cA;
    c_list(i) = c;
%     h_Nlist(i) = N;
    
%     c_temp = 0;
%     for N = N_list
%         fprintf('N = %3d\n', N);
%         case1.N = N;
%         c = case1.solvers(sol_var2{:});
%         if abs(c-c_temp)<1e-8
%             break;
%         elseif ~isnan(c_temp)
%             c_temp = c;
%         end
%         if ~isnan(case1.zc)
%             sol_var2{6}.zL1 = -case1.zc;
%         end
%     end
    [c,~,cA] = case1.solvers(sol_var2{:});
    cA2(i) = cA;
    c_list2(i) = c;
%     h_Nlist2(i) = N;
    
%     c_temp = 0;
%     for N = N_list
%         fprintf('N = %3d\n', N);
%         case1.N = N;
%         c = case1.solvers(sol_var3{:});
%         if abs(c-c_temp)<1e-8
%             break;
%         elseif ~isnan(c_temp)
%             c_temp = c;
%         end
%         if ~isnan(case1.zc)
%             sol_var3{6}.zL1 = -case1.zc;
%         end
%     end
    [c,~,cA] = case1.solvers(sol_var3{:});
    cA3(i) = cA;
    c_list3(i) = c;
%     h_Nlist3(i) = N;
end
toc;

% %% Plot N vs h
% figure;
% plot(h_list,h_Nlist,'-bo');
% hold on;
% plot(h_list,h_Nlist2,'-ro');
% plot(h_list,h_Nlist3,'-go');
% hold off;
% xlabel('$n\lambda$');
% ylabel('$N$');
% legend('eig(A\textbackslash B)','eig(B\textbackslash A)','qz(A,B)','location','northwest');
% grid on;

%% Plot c_i vs h
figure;
plot(h_list,imag(c_list),'-bo');
hold on;
plot(h_list,imag(c_list2),'-ro');
plot(h_list,imag(c_list3),'-go');
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

%% Plot N vs condition number of A (balanced)
figure;
semilogy(h_list,cA1,'-bo');
% hold on;
% semilogy(h_list,cA2,'-ro');
% hold off;
% xlabel('$n\lambda$');
ylabel('cond(A)');
% legend('eig(A\textbackslash B)','qz(A,B)','location','northwest');
grid on;