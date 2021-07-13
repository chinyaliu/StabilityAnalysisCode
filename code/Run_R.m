%% Rayleigh equation solver by eq. finite difference (linear spacing)
%% use derivative matrix (more flexible)
% clear;
%% Initialize
[Test_type,k,N,klist,Nlist,w_FD_eq,count,Testlist] = ToInitialize;
Fr2 = 2.25;
oall = [];
figure;hold on;
for iTest = Testlist
%% Velocity profile
switch Test_type
    case 'N'
        N = iTest;
    case 'k'
        k = iTest;
end
h = 2*pi/real(k);
z = linspace(0,-h,N+1).';
dz = h/N;
a = 0.9988;
sigma = 0.8814;
% zIm = 1i*1/2*sin(pi*(0:N)/N)';
zIm=0;
U = 1 - a*(cosh(sigma*(z+zIm))).^-2;
dU = 2*a*sigma*tanh(sigma*(z+zIm)).*sech(sigma*(z+zIm)).^2;
ddU = 2*sigma^2*a*sech(sigma.*(z+zIm)).^2.*(sech(sigma*(z+zIm)).^2-2*tanh(sigma*(z+zIm)).^2);
[D1,D2] = deri_mat_eqFD(N,dz);
% zbig = ((1:-1:-N-1)/N*h).';
% [D1,D2] = deri_mat_neqFD(N,zbig);
%% Create A (new version)
% N+1 grid point, 2 artificial grid, 1 wave height => "N+4 unknowns"
%%% Kinematic boundary condition
A1_new = [k*U(1), 0, k,zeros(1,N+1)];
%%% Dynamic boundary condition
A2_new = [k/Fr2, 0, -k*dU(1), zeros(1,N+1)] + k*U(1)*D1(1,:);
%%% Governing equation
tmp = ones(N+1,1);
A3_0 = spdiags(tmp,2,N+1,N+4);
A3_new = ((-k^3.*U-k.*ddU).*A3_0).*dz^2 + (k.*U.*D2)*dz^2;
%%% Free slip boundary condition
A4_new = zeros(1,N+4); A4_new(N+3) = 1;
A_new = [A1_new;A2_new;A3_new;A4_new];
%% Create B (new version)
% N+1 grid point, 2 artificial grid, 1 wave height => "N+4 unknowns"
%%% Kinematic boundary condition
B1_new = [1,zeros(1,N+3)];
%%% Dynamic boundary condition
B2_new = D1(1,:);
%%% Governing equation
tmp = ones(N+1,1);
B3_0 = spdiags(tmp,2,N+1,N+4);
B3_new = ((-k^2).*B3_0 + D2)*dz^2;
%%% Free slip boundary condition
B4_new = zeros(1,N+4);
B_new = [B1_new;B2_new;B3_new;B4_new];
%% Generalized eigenvalue problem
fullA_new = full(A_new);
fullB_new = full(B_new);
omega = eig(fullA_new\fullB_new);
omega = 1./omega;
c = omega/k;
c_chosen = filt(c);
if ~isempty(c_chosen)
    oall = [oall;c_chosen*k];
end

% omega_filter = omega(imag(omega) ~= 0);
% [~,ind] = sort(imag(omega_filter),'descend');
% omega_imag_sort = omega_filter(ind);
% count = count + 1;
% if ~isempty(omega_imag_sort)
%     w_FD_eq(count) = omega_imag_sort(1);
% end
disp(['# of grid:',sprintf('% 5d',N),sprintf(', k = %.2f',k)]);
cc = c(abs(c)<10);
gg(1) = scatter(real(cc*k),imag(cc*k),'b');
gg(2) = scatter(real(c_chosen*k),imag(c_chosen*k),'r','filled');
pause(1);
delete(gg);
end 
%%
% ToPlot(Test_type, w_FD_eq, k, Fr2, Nlist, klist,'o-');
% ToSave(Test_type, w_FD_eq, N, k, Fr2, Nlist, klist)

function [D1,D2] = deri_mat_neqFD(N,zbig)
% N => # of derivative points (# of grid points - 2)
% dx => grid size
diffz = abs(diff(zbig));
tmp1 = zeros(N+1,1);
tmp2 = zeros(N+1,1);
tmp3 = zeros(N+1,1);
for i = 1:numel(diffz)-1
    tmp1(i) = diffz(i+1)/(diffz(i)*(diffz(i)+diffz(i+1)));
    tmp2(i) = (diffz(i)-diffz(i+1))/(diffz(i+1)*(diffz(i)));
    tmp3(i) = -diffz(i)/(diffz(i+1)*(diffz(i)+diffz(i+1)));
end
Aplus = spdiags(tmp1,1,N+1,N+4);
A0 = spdiags(tmp2,2,N+1,N+4);
Aminus = spdiags(tmp3,3,N+1,N+4);
D1 = Aplus + A0 + Aminus;


for i = 1:numel(diffz)-1
    tmp1(i) = 2/(diffz(i)*(diffz(i)+diffz(i+1)));
    tmp2(i) = -2/(diffz(i+1)*(diffz(i)));
    tmp3(i) = 2/(diffz(i+1)*(diffz(i)+diffz(i+1)));
end
Bplus = spdiags(tmp1,1,N+1,N+4);
B0 = spdiags(tmp2,2,N+1,N+4);
Bminus = spdiags(tmp3,3,N+1,N+4);
D2 =  Bplus + B0 + Bminus;
end



function [D1,D2] = deri_mat_eqFD(N,dx)
%% Central difference
% N => # of derivative points (# of grid points - 2)
% dx => grid size
tmp = ones(N+1,1);
Aplus = spdiags(tmp,1,N+1,N+4);
A0 = spdiags(tmp,2,N+1,N+4);
Aminus = spdiags(tmp,3,N+1,N+4);
D1 = (Aplus - Aminus)/(2*dx);
D2 = (Aplus - 2*A0 + Aminus)/(dx^2);
%%%%%%%%%%%%%%%%% Test code %%%%%%%%%%%%%%%%%
% N = 16;
% Lx = 2*pi;
% x = linspace(0,Lx,N+1).';
% dx = Lx/N;
% fx = sin(x);
% edfx = -sin(x);
% [D1,D2] = deri_mat_eqFD(N-2,dx);
% D1(:,1) = []; D2(:,1) = [];
% dfx = D2 * flip(fx);
% dfx = flip(dfx);
% plot(x(2:end-1),dfx,x,edfx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function c_chosen = filt(c)
a = 1:length(c); 
crange = ((real(c)-0.0012>-1e-5) & (real(c)-1<=1e-5));
dis = ((real(c)-1).^2 +imag(c).^2)>1e-4;
aa = a(crange&dis);
abch = isoutlier(imag(c(aa)),'movmedian',50);
cgood = abs(c)<10;
aa = [a(dis&~crange&cgood) aa(abch)];
c_chosen = c(aa);
c_chosen = c_chosen(~isinf(c_chosen));
end

function c_chosen = filt2(c)
c_chosen = c((abs(imag(c))>1e-4) & (abs(c)<10));
end