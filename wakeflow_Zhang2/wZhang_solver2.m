function [o, an, cA, errGEP, dob, z, phi, zc, zL] = wZhang_solver2(N,k,h,Re,Fr2,method,alg,balance,zL1,iter)
%% Order of governing equation
switch lower(method(1))
    case {'ray','ray_match'}
        ord = 2;
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
    case 'trefethen'
        [Du,zeta]=cheb(N);
        D(:,:,1) = eye(N+1);
        for i = 2:ord+1
            D(:,:,i) = D(:,:,i-1)*Du;
        end
    otherwise
        error('Invaid method(2) name');
end
BC{1} = permute(D(1,:,:),[3 2 1]);
BC{2} = permute(D(end,:,:),[3 2 1]);
%% Critical height
% syms zz
% g(zz) = finverse(1 - 0.9988*(cosh(0.8814*zz)).^-2);
g = @(zz) (5000*acosh((-2497/(625*(zz - 1)))^(1/2)/2))/4407;
%% Iterate for domain height z_L
zL = zL1;
flag = 0;
for i = 1:11
    %% Differential matrix & Base flow velocity
    w1 = (2/zL).^(0:1:ord);
    w2 = (2/(h-zL)).^(0:1:ord);
    D1 = reshape((reshape(D,[],ord+1).*w1),size(D,1),[],ord+1);
    D2 = reshape((reshape(D,[],ord+1).*w2),size(D,1),[],ord+1);
    [Ubase,z] = baseflow_zhang2(zeta,h,zL);
    BC{3} = Ubase(1,:);
    %% Construct matrix A B
    switch lower(method(1))
        case 'schimd'
            Dall = {D1(3:end-2,:,:), D2(3:end-2,:,:)};
            Ubase = {Ubase(3:N-1,:); Ubase(N+4:end-2,:)};
            [A, B] = matAB_zhang2('d4', N, Dall, Ubase, BC, k, Re, Fr2, w1, w2);
        case 'ray_match'
            Dall = {D1(2:end-1,:,:), D2(1:end-1,:,:)};
            Ubase = {Ubase(2:size(Ubase,1)/2-1,:), Ubase(size(Ubase,1)/2+1:end-1,:)};
            [A, B] = matAB_zhang2(method(1), N, Dall, Ubase, BC, k, Re, Fr2, w1, w2);
        otherwise
            Dall = {D1(2:end-1,:,:), D2(2:end-1,:,:)};
            Ubase = {Ubase(2:size(Ubase,1)/2-1,:), Ubase(size(Ubase,1)/2+2:end-1,:)};
            [A, B] = matAB_zhang2(method(1), N, Dall, Ubase, BC, k, Re, Fr2, w1, w2);
    end
    %% Find eigenvalue(s)
    if strcmp(balance,'y')
        [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
    else
        [o,an,errGEP,cA] = solveGEP(A,B,'max',alg);
        dob = 0;
    end
    ztemp = -double(g(real(o(1))/k));
    if (strcmp(iter,'n')) % known critical height
        break;
    elseif(abs(zL+ztemp) < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', zL);
        break;
    elseif(imag(o(1)) > 0 && i ~= 10) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, zL);
        zL = -ztemp;
    elseif (flag == 1) % didn't converge or growth rate !> 0
        zL = zL1;
        break;
    else % try inflection point z = -0.74708299
        zL = 0.74708299;
        flag = 1;
        fprintf('iter %2d, try inflection pt\n', i);
    end
end
zc = -double(g(real(o)/k));
M = size(D1,2);
phi = [reshape(reshape(permute(D1(:,:,1:3),[1 3 2]),[],M)*an(1:M),[],3);...
    reshape(reshape(permute(D2(:,:,1:3),[1 3 2]),[],M)*an(M+1:end-1),[],3)];
end