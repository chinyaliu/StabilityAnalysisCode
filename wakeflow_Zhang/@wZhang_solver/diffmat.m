function diffmat(obj)
%% Order of governing equation
switch obj.method(1)
    case 'ray'
        obj.ord = 2;
    case 'd4'
        obj.ord = 4;
        if strcmpi(obj.method(3),'d4')
            if strcmpi(obj.method(2),'trefethen')
                    error('Trefethen''s differential method can''t be used for M = N-2\n');
            end
        end
    otherwise
        error('Invalid method(1) name');
end
%% Differential matrix & BC
switch obj.method(2)
    case 'schimd'
        [obj.zeta,obj.Din] = Dcheb(obj.N,obj.ord,obj.method(1));
    case 'trefethen'
        [Du,obj.zeta]=cheb(obj.N);
        obj.Din(:,:,1) = eye(obj.N+1);
        for i = 2:obj.ord+1
            obj.Din(:,:,i) = obj.Din(:,:,i-1)*Du;
        end
    otherwise
        error('Invaid method(2) name');
end
obj.BC{1} = permute(obj.Din(1,:,:),[3 2 1]);
obj.BC{2} = permute(obj.Din(end,:,:),[3 2 1]);
end