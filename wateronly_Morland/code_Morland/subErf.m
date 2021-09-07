classdef subErf < subRay
    methods(Static) 
        function Ur = baseflow(z,ud,delta)
            % error function profile
            Ur(:,1) = ud*erfc(-2*z/delta/sqrt(pi));
            Ur(:,2) = 4*ud*exp(-4*z.^2/pi/delta^2)/pi/delta;
            Ur(:,3) = Ur(:,2).*(-8*z/pi/delta^2);
        end
    end
end