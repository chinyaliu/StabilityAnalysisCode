function [cA,cB,cAb,cBb] = calcond(obj, des, bal, funcN, addvar)
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
end
% Condition number without balancing
cA = cond(A);
cB = cond(B);
% Balancing
if strcmpi(bal,'y')
    n = length(A);
    % balance
    if isreal(A) && isreal(B) && isa(A,'double') && isa(B,'double')
        [Ab,Bb,~,info] = lapack('DGGBAL(h,i,D,i,D,i,i,i,d,D,d,I)',...
            'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
        if info ~= 0
            error('error in dggbal.')
        end
    else
        [Ab,Bb,~,info] = lapack('ZGGBAL(h,i,Z,i,Z,i,i,i,d,D,d,I)',...
            'S',n,A,n,B,n,0,0,zeros(n,1),zeros(n,1),zeros(6*n,1),0);
        if info ~= 0
            error('error in zggbal.')
        end
    end
    cAb = cond(Ab);
    cBb = cond(Bb);
else
    cAb = cA;
    cBb = cB;
end
end      