function [w,P,C] = neoHookeanOLD(F)
%Declare Constants that will be used in neohookean eqns
    lambda_o = 3;
    mu_o = 5;
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
%Calculate the Lagrangian Modulii
    C = zeros(3,3,3,3);
    D = eye(3); %Kronecker Delta
    for k=1:3
    for L=1:3
    C(:,:,k,L) = lambda_o*Finv'*Finv(L,k) + mu_o*D(:,k)*D(:,L)'...
        -(lambda_o*lnJ - mu_o)*Finv(L,:)'*Finv(:,k)';
    end
    end  
end