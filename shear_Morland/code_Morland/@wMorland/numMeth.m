function numMeth(obj,meth)
     obj.setprop('Re',obj.Re);
%     % inviscid
%     obj.method = ["Ray"];
%     obj.ord = 2;
%     obj.subDclass = @subRay;
%     obj.subD = obj.subDclass(); % initialize
    
    obj.method(2) = lower(meth);
    if strcmpi(obj.method(2),'trefethen')
        error('Trefethen''s differential method can''t be used for M = N-2\n');
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

    function [z, D] = shD(obj2)
        if strcmpi(obj.method(1),'d4')
            [z, D] = Dchebnew(obj2.N,obj2.N-2,obj.ord);
        else
            [z, D] = Dchebnew(obj2.N,obj2.N,obj.ord);
        end
    end
    function [z, D] = trD(obj2)
        [Du,z]=cheb(obj2.N);
        D(:,:,1) = eye(obj2.N+1);
        for i = 2:obj.ord+1
            D(:,:,i) = D(:,:,i-1)*Du;
        end
    end
end