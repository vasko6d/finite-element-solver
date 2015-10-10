function [tri,xyz,DOF] = quad_mesh(tri,xyz,DOF)
%Add midpoints to tris and wxy and DOF
    elems = length(tri(:,1));
    L = length(xyz(:,1))+1;
    new_xyz = zeros(3*elems,3);
    new_tri = zeros(elems,3);
    num = 1;
    del = 0;
    for i=1:elems
        n1 = xyz(tri(i,1),:);
        n2 = xyz(tri(i,2),:);
        n3 = xyz(tri(i,3),:);
        n4 = (n1+n2)*.5;
        n5 = (n3+n2)*.5;
        n6 = (n1+n3)*.5;
    %determine if teh point is already in the matrix       
        if i==1
            new_xyz(num,:)   = n4;
            new_xyz(num+1,:) = n5;
            new_xyz(num+2,:) = n6;
            new_tri(i,:) = [L,L+1,L+2];
            del = del+3;
            num = num+3;
     %if it is already stored find its index
        else  
            C = new_xyz(1:num,:);
            ind4 = find(C(:,1)==n4(1)&C(:,2)==n4(2)&C(:,3)==n4(3));
            ind5 = find(C(:,1)==n5(1)&C(:,2)==n5(2)&C(:,3)==n5(3));
            ind6 = find(C(:,1)==n6(1)&C(:,2)==n6(2)&C(:,3)==n6(3));
            if isempty(ind4) 
                new_xyz(num,:) = n4;
                new_tri(i,1) = L+del;
                del = del +1; num = num+1;
            else
                new_tri(i,1) = L+ind4-1;
            end
            if isempty(ind5) 
                new_xyz(num,:) = n5;
                new_tri(i,2) = L+del;
                del = del +1; num = num+1;
            else
                new_tri(i,2) = L+ind5-1;
            end
            if isempty(ind6) 
                new_xyz(num,:) = n6;
                new_tri(i,3) = L+del;
                del = del +1; num = num+1;
            else
                new_tri(i,3) = L+ind6-1;
            end
        end
    end
%cut xyz down to teh correct size
    new_xyz = new_xyz(1:num-1,:);
%update degree of freedoms
    TOL = 10^-8;
    new_DOF = zeros(size(new_xyz));
    dofs = max(max(DOF))+1;
    for i=1:length(new_xyz(:,1))
        if new_xyz(i,3) < TOL
            plot3(new_xyz(i,1),new_xyz(i,2),new_xyz(i,3),'b.');
        else
            new_DOF(i,3) = dofs;
            dofs = dofs+1;
        end
        if new_xyz(i,1) < TOL
            plot3(new_xyz(i,1),new_xyz(i,2),new_xyz(i,3),'r.');
        else
            new_DOF(i,1) = dofs;
            dofs = dofs+1;
        end
        if new_xyz(i,2) < TOL
            plot3(new_xyz(i,1),new_xyz(i,2),new_xyz(i,3),'g.');
        else
            new_DOF(i,2) = dofs;
            dofs = dofs+1;
        end
    end
%Update tri xyz and DOF
    tri = horzcat(tri,new_tri);
    xyz = vertcat(xyz,new_xyz);
    DOF = vertcat(DOF,new_DOF);
%Quadratif mesh has more nodes and DOF
    fprintf('Quadratic Mesh:\n\tNodes:%d\n\tD.O.F:%d\n\n'...
        ,length(xyz(:,1)),dofs-1);
%plot xs on teh middle nodes      
    plot3(new_xyz(:,1),new_xyz(:,2),new_xyz(:,3),'kx');
    %new_tri = new_tri-(L-1)*ones(size(new_tri));
    %trimesh(new_tri(1:elems,:),new_xyz(:,1),new_xyz(:,2),new_xyz(:,3),...
    %    'FaceAlpha',0,'EdgeColor','k');
end