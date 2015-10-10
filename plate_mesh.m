function [tri,xyz,DOF] = plate_mesh(L,W,n_strips,elem_size)
%Determine paramaters for discretizing plate
    L = L/2; W = W/2;
    dL = L/n_strips;
    n = round(W/elem_size);
    dW = W/n;
%Set up approximate size of the structures we will store the data in
    num_nodes = (n_strips+1)*(n+1);
    num_elems = n*2*n_strips;
    xyz = zeros(num_nodes,3);
    tri = zeros(num_elems,3);
%loop through and make nodes
    i=1;
    for y = 0:dL:L
        for x = 0:dW:W
            xyz(i,:) = [x y 0];
            i = i+1;
        end
    end
    
    j = 0;
    k = 0;
    for i = 2:2:num_elems
        tri(i-1,:)  = [1+j+k,2+j+k,2+j+n+k];
        tri(i,:)= [2+j+k,2+j+n+k,3+j+n+k];
        j = j+1;
        if mod(j,n) == 0, k = k+1; end
    end

    
 %Plot all four copies of this eight sphere   
    figure;
    hold on;
    for i=-1:2:1
        for j=-1:2:1
            if i ==1 && j==1 
                trimesh(tri,i*xyz(:,1),j*xyz(:,2),...
                    xyz(:,3),'FaceAlpha',1,'EdgeColor','k');
            else
                    trimesh(tri(1:num_elems,:),i*xyz(:,1),j*xyz(:,2),...
                    xyz(:,3),'FaceAlpha',0,'EdgeColor','b');
            end           
        end
    end
%Do degree of freedom and plot constraints
    DOF = zeros(size(xyz));
    dofs = 1;
    for i=1:length(xyz(:,1))
        if xyz(i,3) == 0
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'b.');
        else
            DOF(i,3) = dofs;
            dofs = dofs+1;
        end
        if xyz(i,1) == 0
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'r.');
        else
            DOF(i,1) = dofs;
            dofs = dofs+1;
        end
        if xyz(i,2) == 0
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'g.');
        else
            DOF(i,2) = dofs;
            dofs = dofs+1;
        end

    end
    daspect([1 1 1]);
    fprintf('1/4 Plate:\n\tNodes:%d\n\tElements:%d\n\tD.O.F:%d\n'...
        ,num_nodes,num_elems,dofs-1);
end

