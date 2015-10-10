function [P1, P2] = plate_iso_def(streach,gauss_pts,elem_typ,mesh_plot)
if nargin<4
    mesh_plot = true;
end
%Make Mesh
    [TRI,XYZ,DOF] = plate_mesh(1,1,2,.5);
    num_elems = length(TRI(:,1));
    if ~mesh_plot
        close;
    end
%figure out our element type
    switch elem_typ
        case 'lin',   shapefunc = @linTri_Na;
        case 'quad',  
            shapefunc = @quadTri_Na;
            [TRI,XYZ,DOF] = quad_mesh(TRI,XYZ,DOF);
        otherwise,    disp('lin or quad elements');
    end
%FUNCTION Equibaxial
    REF = XYZ;
    def = XYZ; 
    DOF(def(:,1) == 0.5 | def(:,2) == 0.5  ,:) = 0;
    def(def(:,1) == 0.5 | def(:,2) == 0.5,:)...
        = streach*def(def(:,1) == 0.5 | def(:,2) == 0.5,:);
        
    DOF = DOF';%-------------IMPORTANT CONVENTION
    num_nodes = length(DOF(1,:));
    residual = 1;
    steps = 0; 
    TOL = 10^-10;
%Find closest stable mesh
    while residual > TOL
        DOF
        [~,Fint,Fext,K] = assemble(REF',def',TRI,...
            [0;0;0],gauss_pts,elem_typ);
        r = Fint-Fext;
        [K, r] = clipper(K,r,DOF);
        residual = sum(abs(r));
        if residual > TOL;
            %solve for update
                u = K\-r;
            %make update
                def = def';
                o = 1; 
                for i=1:num_nodes*3                  
                    if DOF(i) ~= 0;
                        def(i) = def(i) + u(o);
                        o=o+1;
                    end
                end
                def = def';
        end
        steps = steps + 1;
    end
%plot deformed part
if mesh_plot
    for i=-1:2:1
        for j=-1:2:1
            if i ==1 && j==1 
%                 trimesh(TRI(:,1:3),i*def(:,1),j*def(:,2),...
%                     def(:,3),'FaceAlpha',0,'EdgeColor','g');
            else
                          trimesh(TRI(:,1:3),i*def(:,1),j*def(:,2),...
                    def(:,3),'FaceAlpha',0,'EdgeColor','r');
            end     
        end
    end
end
%find 
    for i=1:num_elems
        %fprintf('elem %d:\n',i);
        T = TRI(i,:); %the map of current triangle nodes
        ref_nodes = zeros(3,length(T));
        def_nodes = zeros(3,length(T));
        for a=1:length(T)
            ref_nodes(:,a) = REF(T(a),:)';
            def_nodes(:,a) = def(T(a),:)';
        end
        PandW = quadra_rule(gauss_pts);
        r = PandW(1,1);
        s = PandW(1,2);
        [Na,Na_a] = shapefunc(r,s);
        [G_i, g_i] = get_tans(ref_nodes,def_nodes,Na,Na_a);
        [g_i,~,Tij,~] = enforcePS(G_i,g_i);
        [w,~,~,P1,F] = neoHookean(G_i, g_i);
%         lambda = norm(g_i(:,3));
%         w
%         F
%         P
%         Tij
     %-----------------------------
 if mesh_plot
       for r=0:0.01:1
            s=0;
            [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
            s=1-r;
            [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
        end
        for s=0:0.01:1
            r = 0;
            [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
            plot3(xq,yq,zq,'.r');
        end
 end
     %-----------------------------
    end
        F = [streach 0; 0 streach];
        [w,P2,~,lambda] = PSNeoHookean(F);
%         w
%         lambda
%         P
end