function c = solverE(obj, alg, des, bal, funcN, addvar, errE, loopnum)
    % Set initial critical height guess
    if isfield(addvar,'zL1')
        obj.zc = addvar.zL1;
    else
        obj.zc = 0;
    end
    % Set subdomain positions and number of colloaction points
    obj.setsubd('init', funcN, addvar);
    % Construct matrix A B
    [A,B] = obj.makeAB();
    % De-singularize
    if strcmpi(des,'y')
        er = -200;
        sinB = all(B<eps,2);
        B(sinB,:) = A(sinB,:)/er;
    else
        er = 0;
    end
    % Choose whether to solve with balancing
    n = length(A);
    if strcmpi(bal,'y')
        % balance
        if isreal(A) && isreal(B) && isa(A,'double') && isa(B,'double')
            [A,B,~,info] = lapack('DGGBAL(h,i,D,i,D,i,i,i,d,D,d,I)',...
                'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
            if info ~= 0
                error('error in dggbal.')
            end
        else
            [A,B,~,info] = lapack('ZGGBAL(h,i,Z,i,Z,i,i,i,d,D,d,I)',...
                'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
            if info ~= 0
                error('error in zggbal.')
            end
        end
    end
    c = [];
    for i = 1:loopnum
        op = solEP;
        cp = op(real(op)>0.5*er)/obj.k;
        c = [c; cp(imag(cp)>-5 & imag(cp)<1 & real(cp)<1 & real(cp)>-0.5)];
    end

    function E = randE
        e = randn(n) + 1i*randn(n);
        E = e*errE/norm(e);
    end
    function ev = solEP
        switch alg
            case 'qr'
                M = A\B+randE;
                [~,D] = eig(M);
                ev = 1./diag(D);
            case 'invB'
                M = B\A+randE;
                [~,D] = eig(M);
                ev = diag(D);
        end
    end
end