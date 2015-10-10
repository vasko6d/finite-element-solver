function [W, Fint, Fext, K] = assemble(REF, DEF,TRI,Pext,gauss_pts, elem_typ)
%Super contrived flip the flipped. becuase i changed convention when i
%moved to assembly. Had to add this to use old consistency.
    REF = REF';
    DEF = DEF';
%define basic things needed in loop
    num_nodes = length(REF(:,1));
    num_elems = length(TRI(:,1));
    W = 0;
    Fint = zeros(3,num_nodes);
    Fext = Fint;
    K    = zeros(3,num_nodes,3,num_nodes);
 
 %loop over each element and determine f f and k
    for i=1:num_elems
        %fprintf('elem %d:\n',i);
        T = TRI(i,:); %the map of current triangle nodes
        ref_nodes = zeros(3,length(T));
        def_nodes = zeros(3,length(T));
        for a=1:length(T)
            ref_nodes(:,a) = REF(T(a),:)';
            def_nodes(:,a) = DEF(T(a),:)';
        end
        %Determine if want radia pressure
            if length(Pext) == 1
                [xq,yq,zq]=transform_from_para(def_nodes,[1/3,1/3]...
                    ,elem_typ);
                P = [xq;yq;zq];
                P = P/norm(P);
                P = P*Pext;
            else
                P = Pext;
            end
        [w,fint,fext,k] = WfK(ref_nodes,def_nodes,...
            gauss_pts,elem_typ,P);
    %Update the global arrays,
        W = W +w;
        for a=1:length(T)
            Fint(:,T(a)) = Fint(:,T(a)) + fint(:,a);
            Fext(:,T(a)) = Fext(:,T(a)) + fext(:,a);
            for b=1:length(T)
                K(:,T(a),:,T(b)) = K(:,T(a),:,T(b)) + k(:,a,:,b);
            end           
        end 
    end
    %K = unwrap_K(K);
end