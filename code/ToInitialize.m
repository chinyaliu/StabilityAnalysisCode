function [Test_type,k,N,klist,Nlist,w,count,Testlist,hlist,cellMS...
    ,cellz,cellU,cellV,cellD,wr,L1list]=ToInitialize
%% Initialize
Test_type = 'k';
switch Test_type
    case 'N'
        k = 0.3;
%         k = 5;
%         k = 10;
%         Nlist = (200:200:2000);
        Nlist =[(40:40:180) (200:200:1000)];
%         Nlist = [400,420];
%         Nlist = ones(1,1)*800;
        w = NaN(size(Nlist));
        wr = NaN(size(Nlist));
        L1list = NaN(size(Nlist));
        hlist = NaN(size(Nlist));
        cellMS = cell(size(Nlist));
        cellz = cell(size(Nlist));
        cellU = cell(size(Nlist));
        cellV = cell(size(Nlist));
        cellD = cell(size(Nlist));
        
        Testlist = Nlist;
        klist = []; N = [];
    case 'k'
        N = 1000;
%         klist = linspace(0.01,4,200);
%         klist = 0.01:0.05:0.6;
%         klist = klist - 0.5i;
        klist = 0.1;
        w = NaN(size(klist));
        wr = NaN(size(klist));
        L1list = NaN(size(klist));
        hlist = NaN(size(klist));
        cellMS = cell(size(klist));
        cellz = cell(size(klist));
        cellU = cell(size(klist));
        cellV = cell(size(klist));
        cellD = cell(size(klist));
        Testlist = klist;
        Nlist = []; k = [];
end

count = 0;
end