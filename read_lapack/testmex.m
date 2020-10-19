%{
%%% C
copyfile(fullfile(matlabroot,'extern','examples','refbook','matrixMultiply.c'),'.')
fileattrib('matrixMultiply.c','+w')
mex -v matrixMultiply.c -lmwblas
A = [1 3 5; 2 4 7];
B = [-5 8 11; 3 9 21; 4 0 8];
X = matrixMultiply(A,B);

%%%
copyfile(fullfile(matlabroot,'extern','examples','refbook','matrixDivide.c'),'.')
mex -v matrixDivide.c -lmwlapack
A = [1 2; 3 4];
B = [5; 6];
X = matrixDivide(A,B);

%%%
mex -setup FORTRAN -v
blaslib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft',...
  'libmwblas.lib');
% lapacklib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft',...
%   'libmwlapack.lib');
mex('-v','-largeArrayDims','calldgemm.F',blaslib)
C = calldgemm(A,B);

%%%
mex -largeArrayDims fulltosparse.F loadsparse.F
full = eye(5);
spar = fulltosparse(full)

%}
%%%

N = 4;
A = [1+3i 1+4i 1+5i 1+6i;
    2+2i 4+3i 8+4i 16+5i;
    3+1i 9+2i 27+3i 81+4i;
    4 16+1i 64+2i 256+3i];

B = [1 2+1i 3+2i 4+3i;
    1+1i 4+2i 9+3i 16+4i;
    1+2i 8+3i 27+4i 64+5i;
    1+3i 16+4i 81+5i 256+6i];


[An,Bn,ilo,ihi,lscale,rscale,work,info] = lapack('ZGGBAL(h,i,Z,i,Z,i,I,I,D,D,D,I)',...
    'B',N,A,N,B,N,0,0,zeros(N),zeros(N),zeros(6*N),0);




