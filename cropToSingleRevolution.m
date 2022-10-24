function [Sout] = cropToSingleRevolution(Sin)
    
    Sout = Sin;
    zci = @(v) find(diff(sign(v)));
    ind=zci(Sin{1}(:,2));
    Sout{1} = Sin{1}(1:ind(3),:);
end

