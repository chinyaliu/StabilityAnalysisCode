function [w,xv,dobalance,bakerr,cA] = balanceAB(A,B,modenum,algorithm)
    n = length(A);
    % balance
    if isreal(A) && isreal(B) && isa(A,'double') && isa(B,'double')
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

    cA = cond(A);
    cAb = cond(Ab);
    if cA > cAb
        [w,V] = solveGEP(Ab,Bb,modenum,algorithm);
        % backward transformation to original eigenvector
        xv =  diag(rscale)*V;
        dobalance = 1;
        cA = cAb;
    elseif cA < cAb
        [w,xv] = solveGEP(A,B,modenum,algorithm);
        dobalance = 0;
    end

    if nargout > 3
        bakerr = max(abs(A*xv(:,1)-w(1)*B*xv(:,1)));
    end
end