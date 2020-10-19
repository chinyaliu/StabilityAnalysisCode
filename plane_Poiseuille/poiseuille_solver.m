function [o, an, cA, errGEP, dobalance] = poiseuille_solver(N,k,Re,method,alg,balance)
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
        [z,D] = Dcheb(N,ord,method(1));
    case 'trefethen'
        [Du,z]=cheb(N);
        D(:,:,1) = eye(N+1);
        for i = 2:ord+1
            D(:,:,i) = D(:,:,i-1)*Du;
        end
    otherwise
        error('Invalid method(2) name');
end
BC{1} = permute(D(1,:,:),[3 2 1]);
BC{2} = permute(D(end,:,:),[3 2 1]);
%% Base flow & BC
Ubase = baseflow_poiseuille(z,N);
%% Construct matrix A B
switch lower(method(1))
    case 'schimd'
        [A, B] = matAB_poiseuille('d4', N, D(3:end-2,:,:), Ubase(3:end-2,:,:), BC, k, Re);
    case 'herbert'
        D = [D(1,:,:); D(4:end-3,:,:); D(end,:,:)];
        Ubase = [Ubase(1,:,:); Ubase(4:end-3,:,:); Ubase(end,:,:)];
        [A, B] = matAB_poiseuille('d4', N, D, Ubase, BC, k, Re);
    otherwise
        [A, B] = matAB_poiseuille(method(1), N, D(2:end-1,:,:), Ubase(2:end-1,:,:), BC, k, Re);
end
%% Find eigenvalue(s)
[o,an,dobalance,errGEP,cA] = solveGEP(A,B,1,'y',alg,balance);
end