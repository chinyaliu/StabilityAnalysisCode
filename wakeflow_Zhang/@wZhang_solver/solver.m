function [o, an, cA, errGEP, dob] = solver(obj, zL1, iter, alg, bal)
%% Iterate for domain height z_L
obj.zL = zL1;
flag = 0;
for i = 1:21
      %% Differential matrix & Base flow velocity
    w1 = (2/obj.zL).^(0:1:obj.ord);
    w2 = (2/(obj.h-obj.zL)).^(0:1:obj.ord);
    D1 = reshape((reshape(obj.Din,[],obj.ord+1).*w1),size(obj.Din,1),[],obj.ord+1);
    D2 = reshape((reshape(obj.Din,[],obj.ord+1).*w2),size(obj.Din,1),[],obj.ord+1);
    Ubase = obj.baseflow();
    obj.BC{3} = Ubase(1,:);
      %% Construct matrix A B
    [A, B] = obj.matAB([D1;D2], Ubase);
      %% Find eigenvalue(s)
    if strcmp(bal,'y')
        [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
    else
        if nargout > 2
            [o,an,errGEP,cA] = solveGEP(A,B,1,alg);
            dob = 0;
        else
            [o,an] = solveGEP(A,B,1,alg);
        end
    end
    ztemp = -obj.g(real(o(1))/obj.k);
    if (strcmp(iter,'n')) % known critical height
        break;
    elseif(abs(obj.zL+ztemp) < 1e-9) % converged
        fprintf('converged to zL = %.8f\n', obj.zL);
        break;
    elseif(imag(o(1)) > 0 && i ~= 10) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zL);
        obj.zL = -ztemp;
    elseif (flag == 1 || zL1 == 0.74708299) % didn't converge or growth rate !> 0
        obj.zL = zL1;
        break;
    else % try inflection point z = -0.74708299
        obj.zL = 0.74708299;
        flag = 1;
        fprintf('iter %2d, try inflection pt\n', i);
    end
end
obj.zc = -obj.g(real(o)/obj.k);
if nargout > 1 % calculate modeshape only when eigenvector is needed
    M = size(D1,2);
    obj.z = [0.5*obj.zL*(obj.zeta-1);0.5*(obj.h-obj.zL)*(obj.zeta-1)-obj.zL];
    obj.phi = [reshape(reshape(permute(D1(:,:,1:3),[1 3 2]),[],M)*an(1:M),[],3);...
        reshape(reshape(permute(D2(:,:,1:3),[1 3 2]),[],M)*an(M+1:end-1),[],3)];
end
end      