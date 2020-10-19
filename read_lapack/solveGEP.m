function [w,xv,dobalance,bakerr,cA] = solveGEP(A,B,modenum,includenegative,algorithm,balance)
if strcmp(balance,'y')
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
else
    [w,xv] = QZ(A,B,algorithm,modenum,includenegative);
    bakerr = max(max(abs(A*xv(:,:,1)-w(1)*B*xv(:,:,1))));
    cA = cond(A);
    dobalance = 0;
end
end