classdef subErf < subRay
    methods(Static) 
        function Ur = baseflow(z,varargin)
            % Flow property only has to be input once
            persistent ud delta
            if isempty(ud)||isempty(delta)
                if length(varargin)==2
                    ud = varargin{1};
                    delta = varargin{2};
                else
                    error("Velocity profile property not specified.");
                end
            end
            % error function profile
            Ur(:,1) = ud*erfc(-2*z/delta/sqrt(pi));
            Ur(:,2) = 4*ud*exp(-4*z.^2/pi/delta^2)/pi/delta;
            Ur(:,3) = Ur(:,2).*(-8*z/pi/delta^2);
        end
    end
end