classdef subErf < subRay
    methods 
        function obj = baseflow(obj,ud,delta)
            % error function profile
            Ur(:,1) = ud*erfc(-2*obj.z/delta/sqrt(pi));
            Ur(:,2) = 4*ud*exp(-4*obj.z.^2/pi/delta^2)/pi/delta;
            Ur(:,3) = Ur(:,2).*(-8*obj.z/pi/delta^2);
            obj.U = Ur;
        end
    end
end