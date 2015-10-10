function [forces, displacements] = TL_plate(GS,elem_typ, plother)
%Check args
if nargin <3
    plother = true;
end
%Make the Meshes for the two problems
%    [TRI,XYZ,DOF] = sphere_mesh(1,2,0.75);
    [TRI,REF,DOF] = plate_mesh(.1,.1,1*GS,.05/GS);
    num_elems = length(TRI(:,1));
%figure out our element type
    switch elem_typ
        case 'lin',   
            shapefunc = @linTri_Na;
            gauss_pts = 1;
        case 'quad',  
            shapefunc = @quadTri_Na;
            gauss_pts = 3;
            [TRI,REF,DOF] = quad_mesh(TRI,REF,DOF);
        otherwise,    disp('lin or quad elements');
    end
%Change the z component of the height 
    def = REF;
    L = max(def(:,2));
    def(:,3) = -0.01*cos(def(:,2)*pi/2/L);
%fix degrees of freedom
    DOF = ones(size(DOF));
    DOF(def(:,1) == 0,1) = 0; 
    DOF(def(:,2) == 0,2) = 0;
    DOF(def(:,2) == L,:) = 0;

%main loop
    DOF = DOF';
    num_nodes = length(DOF(1,:));
    
    TOL = 10^-10;
%Find closest stable mesh
steps = 0;
forces = 10:10:10^3;
displacements = zeros(size(forces));
for f = forces;   
        displacements(f/10) = max(abs(def(:,3)));
        steps = 0; 
        residual = 1;
    while residual > TOL
        [~,Fint,Fext,K] = assemble(REF',def',TRI,...
            [0;0;-f],gauss_pts,elem_typ);    
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
    if mod(f, 500) == 0 && plother
    %plot original grid
        plate_mesh(.1,.1,1*GS,.05/GS);
        title('Deformation at One Half Maximum Load')
        if mod(f,1000) == 0
            title('Deformation at Maximum Load')
        end
        grid;
        xlabel('x(m)');ylabel('y(m)');zlabel('z(m)')
    %Get individual triangles. and plot
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
        %Plot it
            if true
                for r=0:0.05:1
                    s=0;
                    [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
                    plot3(xq,yq,zq,'.r');
                    plot3(-xq,yq,zq,'.r');
                    plot3(xq,-yq,zq,'.r');
                    plot3(-xq,-yq,zq,'.r');
                    s=1-r;
                    [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
                    plot3(xq,yq,zq,'.r');
                    plot3(-xq,yq,zq,'.r');
                    plot3(xq,-yq,zq,'.r');
                    plot3(-xq,-yq,zq,'.r');
                end
                for s=0:0.05:1
                    r = 0;
                    [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
                    plot3(xq,yq,zq,'.r');
                    plot3(-xq,yq,zq,'.r');
                    plot3(xq,-yq,zq,'.r');
                    plot3(-xq,-yq,zq,'.r');
                end
            end
        end  
        daspect = [1 1 1];
    end
end
if plother  
    f = figure; hold on;
    title('Transverse Load vs Max Displacement');
    xlabel('Force (N/m^2)'); ylabel('Displacement (m)');
    plot(forces,displacements); grid;
end
end