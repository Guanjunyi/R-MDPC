    function [r_ct] = realcenter(answer,rho)
        cl_ori = unique(answer);
        N_ori = length(cl_ori);
        for i = 1:N_ori
            label = cl_ori(i);
            index = find(answer==label);
            rho_m = max(rho(index));
            index_2 = find(rho(index)==rho_m);
            a_ct = index(index_2(1));
            r_ct(i) = a_ct; 
        end
        
        
        