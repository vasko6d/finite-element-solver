function [K,r] = clipper(K,r,DOF)
    num_nodes = length(DOF(1,:));
    %DOF = reshape(DOF,num_nodes*3,1);
    K = unwrap_K(K);
    r = reshape(r,num_nodes*3,1);
    o = 0;
    for i=1:num_nodes*3
        if DOF(i) == 0;
            K(i-o,:) = []; 
            K(:,i-o) = [];
            r(i-o) = [];
            o = o + 1;
        end
    end
end