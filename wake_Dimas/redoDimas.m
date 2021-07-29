clear all;
if ~contains(path,'code_Dimas;')
    addpath('code_Dimas');
end    
%% Input arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'y';
eigspec = 'all';
N = 1000;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
klist = linspace(0.01,4,400);
ki_list = [0 -0.5 -1 -1.5 -1.7 -2 -2.5 -2.7];
olist = cell(length(ki_list),1);
tic;
for j = 1:length(ki_list)
    ki = klist + ki_list(j)*1i;
    olist{j} = [];
    for i = 1:length(klist)
        k = ki(i);
        fprintf('k = %.2f\n', k);
        [A, B] = makeAB(k, N, h(k), Fr);
        if strcmpi(bal,'y')
            o = balanceAB(A, B, eigspec, alg);
        else
            o = solveGEP(A, B, eigspec, alg);
        end
        o1 = k.*filt(o./k);
        olist{j} = [olist{j}; o1];
    end
end
toc;

%% Plot eigenvalue spectrum
% Read Dimas's results
imagedata = imread('dimas_fr15.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 200;
im2 = cast(im2,'uint8');
imagesc([0,1.4],[0.1,-0.6],im2); % Fr = 1.5
set(gca,'YDir','normal');
hold on;
% Plot results
for i = 1:length(olist)
    o = olist{i};
    textk = sprintf('$k_i=%+1.1f$',ki_list(i));
    plot(real(o),imag(o),'.','Markersize',4,'DisplayName',textk);
end
hold off;
legend('location','southeastoutside');