function [o, an, cA, errGEP, dob] = solver(obj, alg, bal, zL1)
%% Iterate for domain height z_L
obj.zc = zL1;
flag = 0;
for i = 1:21
      %% Differential matrix & Base flow velocity
    w1 = (2/obj.zc).^(0:1:obj.ord);
    w2 = (2/(obj.h-obj.zc)).^(0:1:obj.ord);
    D1 = reshape((reshape(obj.Din,[],obj.ord+1).*w1),size(obj.Din,1),[],obj.ord+1);
    D2 = reshape((reshape(obj.Din,[],obj.ord+1).*w2),size(obj.Din,1),[],obj.ord+1);
    Ubase = obj.baseflow();
      %% Construct matrix A B
    [A, B] = obj.matAB([D1;D2], Ubase);
      %% Find eigenvalue(s)
    if strcmp(bal,'y')
        [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
    else
        if nargout > 2
            [o,an,errGEP,cA] = solveGEP(A,B,'max',alg);
            dob = 0;
        else
            [o,an] = solveGEP(A,B,'max',alg);
        end
    end
    ztemp = -obj.g(real(o(1))/obj.k);
    if(abs(obj.zc+ztemp) < 1e-9) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(imag(o(1)) > 0 && i ~= 20) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        obj.zc = -ztemp;
    else
        fprintf('Didn''t converge.\n');
        break;
    end
end
if nargout > 1 % calculate modeshape only when eigenvector is needed
    M = size(D1,2);
    obj.z = [0.5*obj.zc*(obj.zeta-1);0.5*(obj.h-obj.zc)*(obj.zeta-1)-obj.zc];
    obj.phi = [reshape(reshape(permute(D1(:,:,1:3),[1 3 2]),[],M)*an(1:M),[],3);...
        reshape(reshape(permute(D2(:,:,1:3),[1 3 2]),[],M)*an(M+1:end-1),[],3)];
end
obj.zc = -obj.g(real(o(1))/obj.k);
end      