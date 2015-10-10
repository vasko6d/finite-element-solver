function [] = plot_tans(X,x,point,elem_typ)
%figure out our element type
    switch elem_typ
        case 'lin', shapefunc = @linTri_Na;
        case 'quad',shapefunc = @quadTri_Na;
        otherwise, disp('lin or quad elements');
    end
%Get shape functions
    r = point(1); s=point(2);
    [Na, Na_a] = shapefunc(r,s);
%set up repeted shape function matrixes to make sums easy
    I3 = eye(3);
    N_mat = zeros(3,length(Na)*3);
    Nr_mat = zeros(3,length(Na)*3);
    Ns_mat = zeros(3,length(Na)*3);
    for i=1:length(Na)
        ind = 1+(i-1)*3;
        N_mat(:,ind:ind+2) = I3*Na(i);
        Nr_mat(:,ind:ind+2) = I3*Na_a(i,1);
        Ns_mat(:,ind:ind+2) = I3*Na_a(i,2);
    end  
%Put the the nodal positions in verticle form
    n = length(X(1,:));
    X_line = reshape(X,n*3,1);
    x_line = reshape(x,n*3,1);
%Get the reference tangent vectors
    Xg = N_mat*X_line;
    Gr = Nr_mat*X_line;
    Gs = Ns_mat*X_line;    
    G3 = cross(Gr,Gs); G3 = G3/norm(G3);
%get the deformed tangent vectors
    xg = N_mat*x_line;
    gr = Nr_mat*x_line;
    gs = Ns_mat*x_line;    
    g3 = cross(gr,gs); g3 = g3/norm(g3);
%plot the referece triange with its tangent vectors
    %plot the reference and deforemed only once
    if ~(ishandle(3) && strcmp(get(3, 'type'), 'figure'))            
        figure(3); hold on;
    %plot the reference element
        trimesh([1 2 3],X(1,:),X(2,:),X(3,:),'FaceAlpha',0,'EdgeColor','k');
        plot3(X(1,1),X(2,1),X(3,1),'rx'); %first pt
        plot3(X(1,2),X(2,2),X(3,2),'bx'); %second pt
    %plot the deformed element
        for r=0:0.01:1
            s=0;
            [xq,yq,zq]=transform_from_para(x,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
            s=1-r;
            [xq,yq,zq]=transform_from_para(x,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
        end
        for s=0:0.01:1
            r = 0;
            [xq,yq,zq]=transform_from_para(x,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
        end
    end
%plot the curvilinear coordinates
    plot3(Xg(1),Xg(2),Xg(3),'ko');
    quiver3(Xg(1),Xg(2),Xg(3),Gr(1),Gr(2),Gr(3),'r');
    quiver3(Xg(1),Xg(2),Xg(3),Gs(1),Gs(2),Gs(3),'b');
    quiver3(Xg(1),Xg(2),Xg(3),G3(1),G3(2),G3(3),'k');
    plot3(xg(1),xg(2),xg(3),'ko');
    quiver3(xg(1),xg(2),xg(3),gr(1),gr(2),gr(3),'r');
    quiver3(xg(1),xg(2),xg(3),gs(1),gs(2),gs(3),'b');
    quiver3(xg(1),xg(2),xg(3),g3(1),g3(2),g3(3),'k');
    daspect([1,1,1]);
end