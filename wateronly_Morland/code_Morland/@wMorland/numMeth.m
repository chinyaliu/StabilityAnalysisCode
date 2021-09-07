function numMeth(obj,meth)
    obj.method = ['Ray', lower(meth)];
    obj.ord = 2;
    %% Differential matrix & BC
    switch obj.method(2)
        case 'schimd'
            obj.dm = @shD;
        case 'trefethen'
            obj.dm = @trD;
        otherwise
            error('Invaid method(2) name');
    end

    function [z, D] = shD(obj2)
        [z, D] = Dchebnew(obj2.N,obj2.N,obj.ord);
    end
    function [z, D] = trD(obj2)
        [Du,z]=cheb(obj2.N);
        D(:,:,1) = eye(obj2.N+1);
        for i = 2:obj.ord+1
            D(:,:,i) = D(:,:,i-1)*Du;
        end
    end
end