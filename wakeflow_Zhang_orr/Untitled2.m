%% h = 1L, k = 0.1
% 100
cs(1,1) = -5.69935554321657 + 0.155435654161913i;
cs(1,2) = 0.0943324372553924 - 0.359774643534542i;
cs(1,3) = 7.64805342765973 - 0.0885146263928389i;
% 1000
cs(2,1) = -5.69042132200112 + 0.0157451056762317i;
cs(2,2) = -0.0903635749346843 - 0.0288947409210721i;
cs(2,3) = 7.64522803387940 - 0.00878703878220784i;
% inf
cs(3,1) = -5.69036042186261;
cs(3,2) = 0.0246719832072093 + 0.0344079342992104i;
cs(3,3) = 7.64521877450795;
%% h =6, k = 0.1
% 100
cs(1,1) = -1.35243512853880 + 0.0741596105496452i;
cs(1,2) = 0.0964176159568175 - 0.106243617763662i;
cs(1,3) = 2.92650940472364 - 0.0194903151119640i;
% 1000
cs(2,1) = -1.34599172404507 + 0.00743722414759224i;
cs(2,2) = 0.123010290800953 + 0.0623048508068197i;
cs(2,3) = 2.92617260767161 - 0.00193739708710329i;
% inf
cs(3,1) = -1.34595057448545;
cs(3,2) = 0.137230193360875 + 0.108602511575087i;
cs(3,3) = 2.92617533571192;
%% h = 1L, k = 0.01
% 100
cs(1,1) = -5.69935554321657 + 0.155435654161913i;
cs(1,2) = 0.0943324372553924 - 0.359774643534542i;
cs(1,3) = 7.64805342765973 - 0.0885146263928389i;
% 1000
cs(2,1) = -5.69042132200112 + 0.0157451056762317i;
cs(2,2) = -0.0903635749346843 - 0.0288947409210721i;
cs(2,3) = 7.64522803387940 - 0.00878703878220784i;
% inf
cs(3,1) = -5.69036042186261;
cs(3,2) = 0.0246719832072093 + 0.0344079342992104i;
cs(3,3) = 7.64521877450795;
%% Read Zhang's results & Plot oi vs k
mk = {'o','+','x'};
ln = {'--','-.',':'};
col = get(gca,'colororder');
dn = {'100','1000','$\infty$'};
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,1000,720]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    for j = 1:3
        plot(k,imag(oss{i}(j,:)),'.','markersize',6,'color',col(j,:));
    end
end
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot c real vs imag
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 2
    for j = 1:3
        plot(real(css{i}(j,:)),imag(css{i}(j,:)),'.','Markersize',12,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$c_r$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot real omega
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(oss{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$\omega_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot real c
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(css{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');