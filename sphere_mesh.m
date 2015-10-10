function [tri,xyz,DOF] = sphere_mesh(R,n_strips,elem_size)
%Determine paramaters for discretizing sphere
    elem_size = elem_size;
    dphi = 90/n_strips;
%Set up approximate size of the structures we will store the data in
    approx_num_nodes = ceil(pi/elem_size*n_strips);
    xyz = zeros(approx_num_nodes,3);
    tri = zeros(2*approx_num_nodes,3);
%Begin the first loop that will iterate over various verticle slices of a 
%sphere
    num_elems = 0;
    delta = 0;
    num_nodes = 0;
    for phi=90:-dphi:0
        %Figure out what the radius of this slice is
            r = cosd(phi);
        %The first slice is just a point so treat it seperately
            if phi == 90
                xyz(1,:) = [0,0,1]; prev = [1,1]; num_nodes=num_nodes+1;
        %use a previous strip and current strip approack to have acces to 
        %two rows at a time. Since two will be needed to make triangles
            else
                cur = [sum(prev),0];
            %Decide how elements to actually put i this cut. And the theta 
            %in %degrees between consecutive elements
                num = round(pi*r/2/elem_size);     
                dtheta = 90/num;
            %Iterate through a quarter of a circle
                for theta=0:dtheta:90
                    x = r*cosd(theta);
                    y = r*sind(theta);
                    z = real(sqrt(1-x^2-y^2));
                        if r == 1, z = 0; end
                    a = sum(cur);
                    xyz(a,:) = [x y z];num_nodes=num_nodes+1;
                    cur(2) = cur(2) +1;
                end
            %First Strip is Different sice slice 1 is only one point
                if sum(prev) == 2
                        tri(1,:) = [1, 2, 3];
                        elem = 1;
                        for i=cur(1):sum(cur)-3 
                            tri(i,:) = [1, i+1,i+2]; elem = elem+1;
                        end
                        delta = delta+elem-cur(2);
                        num_elems = num_elems+elem;
                else
             %The general case we ust alternate between nodes to mae 
             %trianges in teh strip. Since cosecutive cuts may have 
             %differet elementshave ot use a compenstation mechanism to 
             %make sure we do notstreach elements. This is what my skip 
             %stsyem is.
                    s = cur(2)/prev(2);
                    skip = mod(s,1);
                    skip_tot = 0;
                    if s >= 2, skip_tot = 0.1; skip = .99; s = 1.90; end               
                    tri(cur(1)-1+delta,:) = [cur(1),prev(1),cur(1)+1];
                        elem = 1;
              %Declare variables used in triangle making. p counts hwat 
              %node we are on in teh previous cut. J is used in conjunction
              %with i to represent what node of the current slice we are 
              %in. k i sneeded to make sure when we skip a top element we 
              %shift numbers so mod 2 will be shifted
                    p = prev(1);
                    j = 0;
                    k = 0;
                    i = cur(1);
              %Actually start main triangle making loop
                    while p < sum(prev)
                        if mod(i-cur(1)+k,2)==0 && skip_tot<1
                            p = p+1;
                            if p == sum(prev)
                                break;
                            end
                            next = p;
                            j = j+1;
                            tri(i+delta,:) = [tri(i-1+delta,2),...
                                tri(i-1+delta,3),next];elem = elem+1;
                        elseif skip_tot >=1 
                                skip_tot = skip_tot-s;
                                k = k+1;
                                next = i+2-j;
                                if next == sum(cur)
                                    break;
                                end                     
                                tri(i+delta,:) = [tri(i-1+delta,3),...
                                    tri(i-1+delta,2),next];elem = elem+1;
                                skip_tot = skip_tot + skip;
                        else
                            next = i+2-j;
                            if next == sum(cur)
                                    break;
                            end
                            tri(i+delta,:) = [tri(i-1+delta,2),...
                                tri(i-1+delta,3),next];elem = elem+1;
                            skip_tot = skip_tot + skip;
                        end
                        i = i+1;
                    end
                    %Delta is needed because we are keeping tract of things 
                    %in the node level adn there are more triangles than  so
                    %nodes need this offset.
                        delta = delta+elem-cur(2);
                        num_elems = num_elems+elem;
                end %sum(prev) =2
                prev = cur;        
            end %End phi = 90
    end
    %num_elems     
    %num_nodes
    tri = tri(1:num_elems,:);
    xyz = R*xyz(1:num_nodes,:);
    %plot3(xyz(:,1),xyz(:,2),xyz(:,3),'o');%plots all nodes with a circle
%Plot all four copies of this eight sphere   
    figure;
    hold on;
    for i=-1:2:1
        for j=-1:2:1
            for k=-1:2:1
                if i ==1 && j==1 && k==1
                    trimesh(tri,i*xyz(:,1),j*xyz(:,2),...
                        k*xyz(:,3),'FaceAlpha',1,'EdgeColor','k');
                else
                    trimesh(tri(1:num_elems,:),i*xyz(:,1),j*xyz(:,2),...
                    k*xyz(:,3),'FaceAlpha',0,'EdgeColor','b');
                end
            end
        end
    end
    TOL = 10^-8;
    DOF = zeros(size(xyz));
    dofs = 1;
    for i=1:length(xyz(:,1)) 
        if xyz(i,1) < TOL
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'r.');
        else
            DOF(i,1) = dofs;
            dofs = dofs+1;
        end
        if xyz(i,2) < TOL
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'g.');
        else
            DOF(i,2) = dofs;
            dofs = dofs+1;
        end
        if xyz(i,3) < TOL
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'b.');
        else
            DOF(i,3) = dofs;
            dofs = dofs+1;
        end
    end
    daspect([1 1 1]);
    fprintf('1/8 Sphere:\n\tNodes:%d\n\tElements:%d\n\tD.O.F:%d\n'...
        ,num_nodes,num_elems,dofs-1);
end