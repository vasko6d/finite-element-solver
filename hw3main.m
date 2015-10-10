close all;
clear all;
clc;
%Make the Meshes for the two problems
%    [TRI,XYZ,DOF] = sphere_mesh(1,2,0.75);
%    [TRI,XYZ,DOF] = plate_mesh(1,1,2,.5);
%plot switches
    II2 = false;
        %mu_o = 5
        %lambda_o = 3;
        %-->changed in neoHookean.m only
    TLconverge = false;
    balloon = true;
    
    single_consistency = false;
        visual_aid = false;
    mesh_consistency = false;
        visual_aid2 = false;
%Determine elemet type ans gauss points 
    elem_typ = 'lin'; 
    %figure out our element type
        switch elem_typ
            case 'lin',   
                shapefunc = @linTri_Na;
                gauss_pts = 1;
                linORquad = 'Linear Element';
            case 'quad',  
                shapefunc = @quadTri_Na;
                gauss_pts = 3;
                %[TRI,XYZ,DOF] = quad_mesh(TRI,XYZ,DOF);
                linORquad = 'Quadratic Element';
            otherwise,    disp('lin or quad elements');
        end
%---------------------CONSISTENCY RAND DEF---------------------------------
if single_consistency
    %use a manual triangle
        X = [-.5 0 0;.5 0 0; 0 1 0];
        deform_typ = 'rand_all';
        [X,x] = gen_deform(X,elem_typ,deform_typ,.1);
        X = X'; x = x';
        %[W,fint,fext,K] = WfK(X,x,gauss_pts,elem_typ,[0;0;0]);
        %K2 = unwrap_K(K);
        if ~visual_aid
            hs = logspace(-12,0,12);
            Es = zeros(2,length(hs));
            for i=1:length(hs)              
                errors = consistency(X,x,gauss_pts,elem_typ,0,hs(i));
                Es(1,i) = errors(1);
                Es(2,i) = errors(2);           
            end
            f = figure; 
            loglog(hs,Es(1,:),hs, hs.^2); hold on;
            title(['Consistency of Internal Force - ',...
                linORquad]);
            xlabel('Perterbation [h]'); ylabel('error');
            saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\'...
                ,'Fcon_s_', elem_typ, '.png']);
            f = figure; 
            loglog(hs,Es(2,:),hs, hs.^2); hold on;
            title(['Consistency of Stiffness - ',...
                linORquad]);
            xlabel('Perterbation [h]'); ylabel('error');
            saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\'...
                ,'Kcon_S_', elem_typ, '.png']);            
    %Visual Aid
        else
            PandW = quadra_rule(gauss_pts); 
            for gp=gauss_pts
                %get what point we are looking at and make a visual aid
                    r = PandW(gp,1);
                    s = PandW(gp,2);
                    plot_tans(X,x,[r,s],elem_typ); %VISUAL AID
            end
        end
end
%-----------------MESH CONSISTENCY---------------------------------
if mesh_consistency
    %Make the Meshes for the two problems
        [TRI,XYZ,DOF] = sphere_mesh(1,2,0.75);
    %    [TRI,XYZ,DOF] = plate_mesh(1,1,2,.5);
    %assemble the stuff with a prescribed displacement
    %     scale = 1.1;
    %     [W,Fint,Fext,K] = assemble(XYZ',scale*XYZ',TRI,...
    %         [0;0;1],gauss_pts,elem_typ);
    %randomize if we chose 
        xyz = XYZ;
        for r=1:30
            c = ceil(rand()*length(XYZ(:,1)));
            e = ceil(rand()*3);
            if DOF(c,e) ~= 0
                m = 2*(rand()-0.5)*0.1;
                xyz(c,e) = xyz(c,e) + m;   
            end  
        end
    %Add midpoints if quadratic
        switch elem_typ
            case 'lin',   
            case 'quad',  
                [TRI,XYZ,DOF] = quad_mesh(TRI,XYZ,DOF);
                [~,xyz,~] = quad_mesh(TRI,xyz,DOF);
            otherwise,    disp('lin or quad elements');
        end
    %consistenct check  
    if ~visual_aid2
            hs = logspace(-12,0,12);
            Es = zeros(2,length(hs));
            for i=1:length(hs)              
                errors = consistency(XYZ',xyz',gauss_pts,elem_typ,TRI,hs(i));
                Es(1,i) = errors(1);
                Es(2,i) = errors(2);           
            end
            f = figure; 
            loglog(hs,Es(1,:),hs, hs.^2); hold on;
            title(['Consistency of Internal Force - ',...
                linORquad]);
            xlabel('Perterbation [h]'); ylabel('error');
            saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\'...
                ,'Fcon_m_', elem_typ, '.png']);
            f = figure; 
            loglog(hs,Es(2,:),hs, hs.^2); hold on;
            title(['Consistency of Stiffness - ',...
                linORquad]);
            xlabel('Perterbation [h]'); ylabel('error');
            saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\'...
                ,'Kcon_m_', elem_typ, '.png'])
    end
    %plot deformed part
    if visual_aid2
        for i=-1:2:1
            for j=-1:2:1
                if i ==1 && j==1 
                    trimesh(TRI(:,1:3),i*xyz(:,1),j*xyz(:,2),...
                        xyz(:,3),'FaceAlpha',0,'EdgeColor','g');
                else
    %                           trimesh(TRI(1:num_elems,:),i*XYZ(:,1),j*XYZ(:,2),...
    %                     XYZ(:,3),'FaceAlpha',0,'EdgeColor','b');
                end     
            end
        end
    end
end
%--------------------------------------------------------------------------
%NOTE
    %K(i,a,k,b) = K(i+3*(a-1),j+3*(b-1))
%------------------------------------------------------------------------
if II2
    maxL = 1;
    [P1,P2,Ls] = plate_iso_def2(maxL,gauss_pts,elem_typ,false);
    f = figure; hold on;
    title('Analytical vs Numerical Solution for Flat Plate Plane Stress');
    xlabel('$\epsilon_{11} = \epsilon_{22}$','Interpreter','Latex'); 
    ylabel('1_s_t P-K stress');
    plot(Ls,P1(1,:),'or');
    plot(Ls,P2(1,:),'b');
    plot(Ls,P1(2,:),'or');
    plot(Ls,P2(2,:),'b');
    legend('Numerical P_{11}','Analytical P_{11}',...
        'Numerical P_{22}','Analytical P_{22}');
    legend('Location', 'NorthEastOutside');
    %saveas(f, 'C:\Users\David\Documents\Latex\261B\hw3\II2_lin.png');
end
%--------------------------------------------------------------------------    
%TL_plate(3,'lin');
if TLconverge
    elem_typ = 'lin';
%Check For Convergence
    if strcmp('lin',elem_typ), 
        gs=[1 2 3 4 5];
    else
        gs=[1 2 3 4];
    end
    Ds = zeros(length(gs),100);
    for i=1:length(gs);
        [Fs,D] = TL_plate(gs(i),elem_typ);
        close all;
        Ds(i,:) = D;  
    end
    f = figure; hold on;
    title('Transverse Load vs Max Displacement');
    xlabel('Force (N/m^2)'); ylabel('Displacement (m)');
    grid; colors = ['k','b','r','g','c'];
    for i =1:length(gs)
        plot(Fs,Ds(i,:),colors(i)); 
    end 
    legend('2 Elements','8 Elements','18 Elements','32 Elements'...
        ,'50 Elements')
    legend('Location','SouthEast');
%     saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\II3_'...
%         ,elem_typ,'_conv1.png']);    
    f = figure; hold on;
    title('Convergence History'); grid;
    xlabel('Elements'); ylabel('Displacement (m)');
    plot(gs.^2*2,Ds(:,end));
%     saveas(f,['C:\Users\David\Documents\Latex\261B\hw3\II3_'...
%         ,elem_typ,'_conv2.png']);
end
%-------------------------------------------------------------------

if balloon
    elem_typ = 'quad';
%Make the Meshes for the two problems
    GS = 3;
    [TRI,REF,DOF] = sphere_mesh(.1,1*GS,pi/2/GS);
    num_elems = length(TRI(:,1));
%figure out our element type
    switch elem_typ
        case 'lin',   shapefunc = @linTri_Na;
            gauss_pts = 1;
        case 'quad',  
            shapefunc = @quadTri_Na;
            gauss_pts = 3;
            [TRI,REF,DOF] = quad_mesh(TRI,REF,DOF);
        otherwise,    disp('lin or quad elements');
    end
% %Change the z component of the height 
    def = REF;
    %Visual Aid
        
%             PandW = quadra_rule(gauss_pts); 
%             for gp=gauss_pts
%                 %get what point we are looking at and make a visual aid
%                     r = PandW(gp,1);
%                     s = PandW(gp,2);
%                     plot_tans(REF',def',[r,s],elem_typ); %VISUAL AID
%             end
%             axis([-1 1 -1 1 -1 1]);
    %randomize
        for r=1:200
            c = ceil(rand()*length(REF(:,1)));
            e = ceil(rand()*3);
            if DOF(c,e) ~= 0
                m = (rand())*0.0001;
                def(c,e) = def(c,e) + m;   
            end  
        end
    
% %main loop
    DOF = DOF';
    num_nodes = length(DOF(1,:)); 
    TOL = 10^-10;
    leaveloop = false;
%Find closest stable mesh
steps = 0;
forces = 10:10:10000;
displacements = zeros(size(forces));
for f = forces; 
        displacements(f/10) = ...
            max(sqrt(def(:,1).^2+def(:,2).^2+def(:,3).^2))/0.1;
        steps = 0; 
        residual = 1;
    while residual > TOL
        [~,Fint,Fext,K] = assemble(REF',def',TRI,...
            f,gauss_pts,elem_typ);    
        r = Fint-Fext;
        [K, r] = clipper(K,r,DOF);
%         rank(K)
%         length(K(:,1))
        if(rank(K) < length(K(:,1))), 
            disp('singular'); 
            leaveloop = true;
            break; 
        end
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
    if mod(f, max(forces)/2) == 0 
    %plot original grid
        sphere_mesh(.1,1*GS,pi/2/GS);
        title('Deformation at One Half Maximum Load')
        if mod(f,max(forces)) == 0
            title('Deformation at Maximum Load')
        end
        grid;
        xlabel('x(m)');ylabel('y(m)');zlabel('zm)')
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
                    plot3(xq,yq,-zq,'.r');
                    plot3(-xq,yq,-zq,'.r');
                    plot3(xq,-yq,-zq,'.r');
                    plot3(-xq,-yq,-zq,'.r');
                    s=1-r;
                    [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
                    plot3(xq,yq,zq,'.r');
                    plot3(-xq,yq,zq,'.r');
                    plot3(xq,-yq,zq,'.r');
                    plot3(-xq,-yq,zq,'.r');
                    plot3(xq,yq,-zq,'.r');
                    plot3(-xq,yq,-zq,'.r');
                    plot3(xq,-yq,-zq,'.r');
                    plot3(-xq,-yq,-zq,'.r');
                end
                for s=0:0.05:1
                    r = 0;
                    [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
                    plot3(xq,yq,zq,'.r');
                    plot3(-xq,yq,zq,'.r');
                    plot3(xq,-yq,zq,'.r');
                    plot3(-xq,-yq,zq,'.r');
                    plot3(xq,yq,-zq,'.r');
                    plot3(-xq,yq,-zq,'.r');
                    plot3(xq,-yq,-zq,'.r');
                    plot3(-xq,-yq,-zq,'.r');
                    
                end
            end
        end  
        daspect = [1 1 1];
    end
    if leaveloop
        break;
    end
end
    
f = figure; hold on;
title('Internal Pressure vs Stretch Ratio');
xlabel('Pressure (N/m^2)'); ylabel('Stretch Ratio r/R');
plot(forces,displacements); grid;
end

















% %Make the Meshes for the two problems
% %    [TRI,XYZ,DOF] = sphere_mesh(1,2,0.75);
%     [TRI,REF,DOF] = plate_mesh(.1,.1,2,.05);
%     num_elems = length(TRI(:,1));
% %figure out our element type
%     switch elem_typ
%         case 'lin',   shapefunc = @linTri_Na;
%         case 'quad',  
%             shapefunc = @quadTri_Na;
%             [TRI,REF,DOF] = quad_mesh(TRI,REF,DOF);
%         otherwise,    disp('lin or quad elements');
%     end
% %Change the z component of the height 
%     def = REF;
%     L = max(def(:,2));
%     def(:,3) = -0.01*cos(def(:,2)*pi/2/L);
% %fix degrees of freedom
%     DOF = ones(size(DOF));
%     DOF(def(:,1) == 0,1) = 0; 
%     DOF(def(:,2) == 0,2) = 0;
%     DOF(def(:,2) == L,:) = 0;
% 
% %main loop
%     DOF = DOF';
%     num_nodes = length(DOF(1,:));
%     
%     TOL = 10^-10;
% %Find closest stable mesh
% steps = 0;
% forces = 10:10:10^3;
% displacements = zeros(size(forces));
% for f = forces;   
%         displacements(f/10) = max(abs(def(:,3)));
%         steps = 0; 
%         residual = 1;
%     while residual > TOL
%         [~,Fint,Fext,K] = assemble(REF',def',TRI,...
%             [0;0;-f],gauss_pts,elem_typ);    
%         r = Fint-Fext;
%         [K, r] = clipper(K,r,DOF);
%         residual = sum(abs(r));
%         if residual > TOL;
%             %solve for update
%                 u = K\-r;
%             %make update
%                 def = def';
%                 o = 1; 
%                 for i=1:num_nodes*3                  
%                     if DOF(i) ~= 0;
%                         def(i) = def(i) + u(o);
%                         o=o+1;
%                     end
%                 end
%                 def = def';
%         end
%         steps = steps + 1;
%     end
%     if mod(f, 500) == 0 
%     %plot original grid
%         plate_mesh(.1,.1,2,.05);
%         title('Deformation at One Half Maximum Load')
%         if mod(f,1000) == 0
%             title('Deformation at Maximum Load')
%         end
%         grid;
%         xlabel('x(m)');ylabel('y(m)');zlabel('zm)')
%     %Get individual triangles. and plot
%         for i=1:num_elems
%             %fprintf('elem %d:\n',i);
%             T = TRI(i,:); %the map of current triangle nodes
%             ref_nodes = zeros(3,length(T));
%             def_nodes = zeros(3,length(T));
%             for a=1:length(T)
%                 ref_nodes(:,a) = REF(T(a),:)';
%                 def_nodes(:,a) = def(T(a),:)';
%             end
%             PandW = quadra_rule(gauss_pts);
%             r = PandW(1,1);
%             s = PandW(1,2);
%             [Na,Na_a] = shapefunc(r,s); 
%         %Plot it
%             if true
%                 for r=0:0.05:1
%                     s=0;
%                     [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%                     plot3(xq,yq,zq,'.r');
%                     plot3(-xq,yq,zq,'.r');
%                     plot3(xq,-yq,zq,'.r');
%                     plot3(-xq,-yq,zq,'.r');
%                     s=1-r;
%                     [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%                     plot3(xq,yq,zq,'.r');
%                     plot3(-xq,yq,zq,'.r');
%                     plot3(xq,-yq,zq,'.r');
%                     plot3(-xq,-yq,zq,'.r');
%                 end
%                 for s=0:0.05:1
%                     r = 0;
%                     [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%                     plot3(xq,yq,zq,'.r');
%                     plot3(-xq,yq,zq,'.r');
%                     plot3(xq,-yq,zq,'.r');
%                     plot3(-xq,-yq,zq,'.r');
%                 end
%             end
%         end  
%         daspect = [1 1 1];
%     end
% end
%     
% f = figure; hold on;
% title('Transverse Load vs Max Displacement');
% xlabel('Force (N/m^2)'); ylabel('Displacement (m)');
% plot(forces,displacements); grid;

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% %FUNCTION Equibaxial-----------------------------------------------------
%     REF = XYZ;
%     def = XYZ;
%     streach = 1.1; 
%     DOF(def(:,1) == 0.5 | def(:,2) == 0.5  ,:) = 0;
%     def(def(:,1) == 0.5 | def(:,2) == 0.5,:)...
%         = streach*def(def(:,1) == 0.5 | def(:,2) == 0.5,:);
%         
%     DOF = DOF';%-------------IMPORTANT CONVENTION
%     num_nodes = length(DOF(1,:));
%     residual = 1;
%     steps = 0; 
%     TOL = 10^-10;
%     while residual > TOL
%         steps
%         [W,Fint,Fext,K] = assemble(REF',def',TRI,...
%             [0;0;0],gauss_pts,elem_typ);    
%         r = Fint-Fext;
%         [K, r] = clipper(K,r,DOF);
%         residual = sum(abs(r));
%         if residual > TOL;
%             %solve for update
%                 u = K\-r;
%             %make update
%                 def = def';
%                 o = 1; 
%                 for i=1:num_nodes*3                  
%                     if DOF(i) ~= 0;
%                         def(i) = def(i) + u(o);
%                         o=o+1;
%                     end
%                 end
%                 def = def';
%         end
%         steps = steps + 1;
%     end
% %plot deformed part
%     num_elems = length(TRI(:,1));
%     for i=-1:2:1
%         for j=-1:2:1
%             if i ==1 && j==1 
% %                 trimesh(TRI(:,1:3),i*def(:,1),j*def(:,2),...
% %                     def(:,3),'FaceAlpha',0,'EdgeColor','g');
%             else
%                           trimesh(TRI(:,1:3),i*def(:,1),j*def(:,2),...
%                     def(:,3),'FaceAlpha',0,'EdgeColor','r');
%             end     
%         end
%     end
% 
%     for i=1:num_elems
%         %fprintf('elem %d:\n',i);
%         T = TRI(i,:); %the map of current triangle nodes
%         ref_nodes = zeros(3,length(T));
%         def_nodes = zeros(3,length(T));
%         for a=1:length(T)
%             ref_nodes(:,a) = REF(T(a),:)';
%             def_nodes(:,a) = def(T(a),:)';
%         end
%         PandW = quadra_rule(gauss_pts);
%         r = PandW(1,1);
%         s = PandW(1,2);
%         [Na,Na_a] = shapefunc(r,s);
%         [G_i, g_i] = get_tans(ref_nodes,def_nodes,Na,Na_a);
%         [~,Tau,~,~,~] = neoHookean(G_i, g_i);
%         [g_i,~,Tij,~] = enforcePS(G_i,g_i);
%         [w,~,~,P,F] = neoHookean(G_i, g_i);
%         lambda = norm(g_i(:,3));
%         w
%         F
%         P
%         Tij
%      %-----------------------------
%        for r=0:0.01:1
%             s=0;
%             [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%             plot3(xq,yq,zq,'.r');
%             s=1-r;
%             [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%             plot3(xq,yq,zq,'.r');
%         end
%         for s=0:0.01:1
%             r = 0;
%             [xq,yq,zq]=transform_from_para(def_nodes,[r,s],elem_typ);
%             plot3(xq,yq,zq,'.r');
%         end
%      %-----------------------------
%     end
%         F = [streach 0; 0 streach];
%         [w,P,~,lambda] = PSNeoHookean(F);
%         w
%         lambda
%         P
 
%----------original assemble----------------
% %Make global structures for fint and fext and kiakb
%     num_nodes = length(XYZ(:,1));
%     num_elems = length(TRI(:,1));
%     W = 0;
%     Fint = zeros(3,num_nodes);
%     Fext = Fint;
%     K    = zeros(3,num_nodes,3,num_nodes);
%     for i=1:num_elems
%         T = TRI(i,:); %the map of current triangle nodes
%         nodes = zeros(3,length(T));
%         for a=1:length(T)
%             nodes(:,a) = XYZ(T(a),:)';
%         end
%         [w,fint,fext,k] = WfK(nodes,scale*nodes,gauss_pts,elem_typ,[0;0;1]);
%         W = W +w;
%         for a=1:length(T)
%             Fint(:,T(a)) = Fint(:,T(a)) + fint(:,a);
%             Fext(:,T(a)) = Fext(:,T(a)) + fext(:,a);
%             for b=1:length(T)
%                 K(:,T(a),:,T(b)) = K(:,T(a),:,T(b)) + k(:,a,:,b);
%             end           
%         end 
%     end
%     K = unwrap_K(K);















