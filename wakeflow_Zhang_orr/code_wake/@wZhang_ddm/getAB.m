function [A,B] = getAB(obj, funcN, addvar)
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
end