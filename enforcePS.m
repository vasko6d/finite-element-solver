function [g_i,w,Tij,Cijkl] = enforcePS(G_i,g_i)
%initial lambda
    lambda = norm(g_i(:,3)); %should be 1 initially
    g_i(:,3) = g_i(:,3)/lambda;
%Declare the Toleranceofour NR
    TOL = 10^-8;
%Some initial conditions for the NR
    g_i2 = g_i;
    [Gi,gi] = get_duals(G_i,g_i2);
    [w,Tau,CIJKL] = neoHookean(G_i,g_i2);
    Tij = get_Tij(Tau,gi);
    Cijkl = get_Cijkl(CIJKL,Gi);
    counter  = 0;
%loop until Tij(3,3) is zero
    while abs(Tij(3,3)) > TOL
        counter = counter +1;   
        %disp(counter);
        d_lambda =  - Tij(3,3)/2/lambda/Cijkl(3,3,3,3);
        lambda = lambda + d_lambda;
        g_i2(:,3) = lambda*g_i(:,3); 
        [Gi,gi] = get_duals(G_i,g_i2);
        [w,Tau,CIJKL] = neoHookean(G_i,g_i2);
        Cijkl = get_Cijkl(CIJKL,Gi);
        Tij = get_Tij(Tau,gi);
        if counter > 10, disp('NR not converging'); end
        if counter > 50, break; end
        %Tij(3,3);
    end
%scale our director correctly
    g_i(:,3) = lambda*g_i(:,3); 
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
%Get the stresses in teh correct reference
    function [Tij] = get_Tij(Tau,gi)
        Tij = zeros(3,3);
        for i=1:3
        for j=1:3
        for I=1:3
        for J=1:3
            Tij(i,j) = Tij(i,j) + gi(I,i)*Tau(I,J)*gi(J,j);
        end
        end        
        end
        end 
    end
%Get tangent moduli in teh correct reference
    function [Cijkl] = get_Cijkl(CIJKL,Gi)
        Cijkl = zeros(3,3,3,3);
%         for i=1:3
%         for j=1:3
%         for K=1:3
%         for L=1:3
%         for k=1:3
%         for l=1:3
%         for J=1:3
%         for I=1:3           
%             Cijkl(i,j,k,l) = Cijkl(i,j,k,l)+CIJKL(I,J,K,L)*Gi(I,i)*Gi(J,j)...
%                 *Gi(K,k)*Gi(L,l);
%         end
%         end
%         end
%         end
%         end
%         end
%         end
%         end
        for k=1:3
        for l=1:3
        for J=1:3
        for I=1:3
            tx = reshape(CIJKL(I,J,:,:),3,3).*(Gi(:,k)*Gi(:,l)');
            T = sum(sum(tx));
            M = (Gi(I,:)'*Gi(J,:));               
            Cijkl(:,:,k,l) = Cijkl(:,:,k,l) + M*T;
        end
        end
        end
        end
    end
end