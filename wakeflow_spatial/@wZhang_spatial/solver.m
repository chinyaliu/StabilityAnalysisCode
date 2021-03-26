function [k, an, cA, errGEP, dob] = solver(obj, alg, bal)
%% Differential matrix & Base flow velocity
Ubase = obj.baseflow();
%% Construct matrix A
[A0, A1, A2, A3] = obj.matAB(obj.Din, Ubase);
A = {A0, A1, A2, A3};
%% Find eigenvalue(s)
if strcmp(bal,'y')
    %%[o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
else
    dob = 0;
    if nargout > 2
        [k,an,errGEP,cA] = solvePoly(A,'all',alg);
    else
        [k,an] = solvePoly(A,'all',alg);
    end
end
% if nargout > 1 % calculate modeshape only when eigenvector is needed
%    M = size(obj.Din,2);
%    obj.z = [0.5*obj.zc*(obj.zeta-1);0.5*(obj.h-obj.zc)*(obj.zeta-1)-obj.zc];
%    obj.phi = [reshape(reshape(permute(D1(:,:,1:3),[1 3 2]),[],M)*an(1:M),[],3);...
%        reshape(reshape(permute(D2(:,:,1:3),[1 3 2]),[],M)*an(M+1:end-1),[],3)];
% end
% obj.zc = -obj.g(obj.o/real(k));
end      