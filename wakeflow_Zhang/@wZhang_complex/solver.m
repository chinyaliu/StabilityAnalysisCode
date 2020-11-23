function [o, an, cA, errGEP, dob] = solver(obj, alg, bal)
%% Differential matrix & Base flow velocity
% w = 1./((0.5*obj.h*(1-2*1i*obj.del*obj.zeta)).^(0:1:obj.ord));
w = 1./((0.5*obj.h*(1-0.5*pi*1i*obj.del*sinpi(obj.zeta/2))).^(0:1:obj.ord));
D = permute(permute(obj.Din,[1 3 2]).*w, [1 3 2]);
Ubase = obj.baseflow();
%% Construct matrix A B
[A, B] = obj.matAB(D, Ubase);
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
% obj.zc = -obj.g(real(o)/obj.k);
if nargout > 1 % calculate modeshape only when eigenvector is needed
    M = size(D,2);
    obj.phi = reshape(reshape(permute(D(:,:,1:3),[1 3 2]),[],M)*an(1:M),[],3);
end
end