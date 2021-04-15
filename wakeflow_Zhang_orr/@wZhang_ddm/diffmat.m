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
        obj.dm = @shD;
    case 'trefethen'
        obj.dm = @trD;
    otherwise
        error('Invaid method(2) name');
end
    function shD(obj2)
        [obj2.zeta, obj2.Din] = Dcheb(obj2.N,obj.ord,obj.method(1));
    end
    function trD(obj2)
        [Du,obj2.zeta]=cheb(obj2.N);
        Din(:,:,1) = eye(obj2.N+1);
        for i = 2:obj.ord+1
            Din(:,:,i) = Din(:,:,i-1)*Du;
        end
        obj2.Din = Din;
    end
end