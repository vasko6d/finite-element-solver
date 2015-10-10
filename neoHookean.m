function [w,Tau,CIJKL,P,F] = neoHookean(G_i,g_i)
%Calculate F
    F = get_F(G_i,g_i);
%Declare Constants that will be used in neohookean eqns
    lambda_o = 40*10^5;
    mu_o = 4*10^5;
%Right Cauchy Green Deformation Tensor
%and other usefull constants
    RCGDT = F.'*F;
    J = det(F);
    I_1 = trace(RCGDT);
    Finv = F^-1;
    lnJ = log(J);
%strain energy density   
    w = lambda_o/2*(lnJ^2) - mu_o*lnJ + mu_o/2*(I_1 - 3);
%1st PK stress
    P = (lambda_o*lnJ-mu_o)*Finv'+mu_o*F;
    Tau = P*F';
%Calculate the Lagrangian Modulii
    CiJkL = zeros(3,3,3,3);
    CIJKL = CiJkL;
    S = Finv*P;
    D = eye(3); %Kronecker Delta
    for k=1:3
    for L=1:3
    CiJkL(:,:,k,L) = lambda_o*Finv'*Finv(L,k) + mu_o*D(:,k)*D(:,L)'...
        -(lambda_o*lnJ - mu_o)*Finv(L,:)'*Finv(:,k)';
    end
    end
%convert it to CIJKL
    for I=1:3
    for J=1:3
    for K=1:3
    for L=1:3
    for k=1:3
    for i=1:3
    CIJKL(I,J,K,L) = CIJKL(I,J,K,L) + 0.5*Finv(I,i)*Finv(K,k)...
        *(CiJkL(i,J,k,L)-D(i,k)*S(J,L));
    end
    end 
    end
    end
    end
    end
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
    function [F] = get_F(G_i,g_i)
    %Back Calculate Lamda and a3
        lambda = norm(g_i(:,3));
        a_3 = g_i(:,3)/lambda;
        g_i(:,3) = a_3;
    %get the metrics
        G_ij = G_i'*G_i;   %g_ij = g_i'*g_i;
        Gij = G_ij^-1;     %gij = g_ij^-1;
    %Dual Vectors
        Gi = G_i*Gij';      %gi = g_i*gij';
    %Deformation Gradient
        F = g_i(:,1)*Gi(:,1)'+g_i(:,2)*Gi(:,2)'+lambda*a_3*Gi(:,3)';
    end
end