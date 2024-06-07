function [] = startParPool(n_cpu)
    flag = 1;
    while(flag)
        try
            p = gcp('nocreate');
            if ~exist('n_cpu','var')
                if isempty(p); parpool('local'); end
            else
                if isempty(p); parpool('local',n_cpu); end
            end
            flag = 0;
        catch
            flag = 1;
        end
    
    end
end