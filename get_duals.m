function [Gi,gi] = get_duals(G_i,g_i)
    %get the metrics
        G_ij = G_i'*G_i;   
        g_ij = g_i'*g_i;
        Gij = G_ij^-1;     
        gij = g_ij^-1;
    %Dual Vectors
        Gi = G_i*Gij';      
        gi = g_i*gij';
end