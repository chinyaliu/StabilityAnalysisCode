function [o, an, cA, errGEP] = wZhang_solver(N,k,h,Re,Fr2,method,alg)
%% Differential order of solving method
switch lower(method(1))
    case 'd2'
        ord = 2;
    case 'uw'
        ord = 3;
    case 'd4'
        ord = 4;
        if strcmpi(method(3),'d4')
            if strcmpi(method{2},'trefethen')
                    error('Trefethen''s differential method can''t be used for M = N-2\n');
            end
        end
        method(1) = method(3);
    otherwise
        error('Invalid method(1) name');
end
%% Differential matrix & BC
switch lower(method(2))
    case 'schimd'
        [zeta,D] = Dcheb(N,ord,method(1));
        D = reshape((reshape(D,[],ord+1).*((2/h).^(0:1:ord))),size(D,1),[],ord+1);
    case 'trefethen'
        [Du,zeta]=cheb(N);
        D(:,:,1) = eye(N+1);
        for i = 2:ord+1
            D(:,:,i) = (2/h)*D(:,:,i-1)*Du;
        end
    otherwise
        error('Invalid method(2) name');
end
BC{1} = permute(D(1,:,:),[3 2 1]);
BC{2} = permute(D(end,:,:),[3 2 1]);
%% Base flow & BC
Ubase = baseflow_zhang(zeta,h);
BC{3} = Ubase(1,:);
%% Construct matrix A B
if isinf(Re)
    [A, B] = matAB_zhang('ray', N, D(2:end-1,:,:), Ubase(2:end-1,:), BC, k, Re, Fr2);
else
    switch lower(method(1))
        case 'schimd'
            [A, B] = matAB_zhang('d4', N, D(3:end-2,:,:), Ubase(3:end-2,:), BC, k, Re, Fr2);
        case 'herbert'
            D = [D(1,:,:); D(4:end-3,:,:); D(end,:,:)];
            Ubase = [Ubase(1,:); Ubase(4:end-3,:); Ubase(end,:)];
            [A, B] = matAB_zhang('d4', N, D, Ubase, BC, k, Re, Fr2);
        otherwise
            [A, B] = matAB_zhang(method(1), N, D(2:end-1,:,:), Ubase(2:end-1,:), BC, k, Re, Fr2);
    end
end
%% Find eigenvalue(s)
[o,an,errGEP,cA] = solveGEP(A,B,'max',alg);
end