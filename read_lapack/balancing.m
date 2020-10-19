function [w,xv,dobalance,bakerr,cA] = balancing(A,B,modenum,includenegative,algorithm)

n = length(A);
% balance
if isdouble(fi(A)) && isdouble(fi(B))
    [Ab,Bb,rscale,info] = lapack('DGGBAL(h,i,D,i,D,i,i,i,d,D,d,I)',...
        'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
    if info ~= 0
        error('error in dggbal.')
    end
else
    [Ab,Bb,rscale,info] = lapack('ZGGBAL(h,i,Z,i,Z,i,i,i,d,D,d,I)',...
        'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
    if info ~= 0
        error('error in zggbal.')
    end
end

% qz algorithm
% [w1,w2,~,V1,V2,~,~,e2] = QZ_filter(Ab,Bb,'eig');
% [alpha,beta,info] = lapack('ZGGEV(h,h,i,z,i,z,i,Z,Z,z,i,z,i,z,i,d,I)',...
%     'V','V',n,Ab,n,Bb,n,zeros(n,1),zeros(n,1),zeros(n),n,zeros(n),n,zeros(2*n,1),2*n,zeros(8*n,1),0);
if cond(A) > cond(Ab)
    cA = cond(Ab);
    [w,V] = QZ(Ab,Bb,algorithm,modenum,includenegative);
    % backward transformation to original eigenvector
    xv = nan(size(V)); % preallocate
    for i = 1:modenum
        xv(:,:,i) = diag(rscale)*V(:,:,i);
    end
    dobalance = 1;
elseif cond(A) < cond(Ab)
    cA = cond(A);
    [w,V] = QZ(A,B,algorithm,modenum,includenegative);
    xv(:,:,1:modenum) = V(:,:,1:modenum);
    dobalance = 0;
end

bakerr = max(max(abs(A*xv(:,:,1)-w(1)*B*xv(:,:,1))));

end