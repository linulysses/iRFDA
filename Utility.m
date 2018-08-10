classdef Utility
    methods (Static)
        
        function b = Contains(s,sarray)
            for k = 1:length(sarray)
                if strcmp(s,sarray{k})
                    b = true;
                    return;
                end
            end
            b = false;
        end
        
    end
    
end