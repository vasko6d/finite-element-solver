function [W, fint,fext, K] = WfK(X, x, gauss_pts,elem_typ,Pext,skip)
%figure out our element type
    switch elem_typ
        case 'lin',   shapefunc = @linTri_Na;
        case 'quad',  shapefunc = @quadTri_Na;
        otherwise,    disp('lin or quad elements');
    end
%flag to save time
    if nargin < 6,skip = false; end
%shape functions
    PandW = quadra_rule(gauss_pts); 
    nn = length(X(1,:));
    H = 0.001;
    mu = 1;
    fint = zeros(3,nn);
    fext = fint;
    K = zeros(3,nn,3,nn);
    W = 0;
%Looop over quadrature points
    for gp=1:gauss_pts
        %get what point we are looking at and make a visual aid
            r = PandW(gp,1);
            s = PandW(gp,2);
            %plot_tans(X,x,[r,s],elem_typ); %VISUAL AID
        %Get the tangent vectors and other useful qunatities
            [Na,Na_a] = shapefunc(r,s);
            [G_i,g_i] = get_tans(X,x,Na,Na_a);
            sqrtA = norm(cross(G_i(:,1),G_i(:,2)));
            [g_i,w,Tij,Cijkl] = enforcePS(G_i,g_i);
        %Get the stress resultant
            na = g_i*Tij(:,1:2);
        %Strain Energy
            W = W + 0.5*PandW(gp,3)*sqrtA*H*mu*w;
        %internal forces
            ftemp = na*Na_a'*mu*H*sqrtA;
            fint = fint + 0.5*PandW(gp,3)*ftemp;
        %exteral forces NOTE: no mu or h here?
            ftemp = Pext*Na*sqrtA;
            fext = fext + 0.5*PandW(gp,3)*ftemp;
        %tangent stiffness
        if ~skip 
        %Get C_2D
            C1 = Cijkl(1:2,1:2,3,3);
            C2 = reshape(Cijkl(3,3,1:2,1:2),2,2)/Cijkl(3,3,3,3);
            Cd = zeros(2,2,2,2);
                Cd(:,:,1,1) = C1*C2(1,1);
                Cd(:,:,2,1) = C1*C2(2,1);
                Cd(:,:,1,2) = C1*C2(1,2);
                Cd(:,:,2,2) = C1*C2(2,2);
            C_2D = Cijkl(1:2,1:2,1:2,1:2) - Cd;
        %actual stiffness
            Kt = zeros(3,nn,3,nn);
            D = eye(3);%kronecker delta
            for i=1:3
            for k=1:3
            for A=1:nn
            for B=1:nn
            for a=1:2
            for b=1:2
                Kt(i,A,k,B) = Kt(i,A,k,B) + Tij(a,b)*D(i,k)...
                    *Na_a(A,a)*Na_a(B,b)*sqrtA;
            for m=1:2
            for n=1:2
                Kt(i,A,k,B) = Kt(i,A,k,B) + 2*C_2D(a,b,m,n)...
                    *g_i(i,b)*g_i(k,n)*Na_a(A,a)*Na_a(B,m)*sqrtA;
            end
            end
            end
            end
            end
            end
            end
            end
            K = K + 0.5*PandW(gp,3)*Kt*H*mu;
        end
    end
end