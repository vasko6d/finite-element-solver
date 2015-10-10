function [w,P,C,lambda] = PSNeoHookean(F)
%make F 3D
    Ffull = zeros(3,3);
    Ffull(1:2,1:2) = F;
    Ffull(3,3) = 1;
%calculare lambda F(3,3)
    d_lambda = 1;
    tol = 10^-10;
    while abs(d_lambda) > tol
        [~,P,C] = neoHookeanOLD(Ffull);
        d_lambda = - P(3,3)/C(3,3,3,3);
        Ffull(3,3) = Ffull(3,3) + d_lambda;
    end
    lambda = Ffull(3,3);
%Declare Constants that will be used in neohookean eqns
    lambda_o = 3;
    mu_o = 5;
%Right Cauchy Green Deformation Tensor
%and other usefull constants
    RCGDT = Ffull.'*Ffull;
    J = det(Ffull);
    I_1 = trace(RCGDT);
    Finv = F(1:2,1:2)^-1;
    lnJ = log(J);
%strain energy density   
    w = lambda_o/2*(lnJ^2) - mu_o*lnJ + mu_o/2*(I_1 - 3);
%1st PK stress
    P = (lambda_o*lnJ-mu_o)*Finv'+mu_o*F;
%Calculate the Lagrangian Modulii
    C = zeros(2,2,2,2);
    S = eye(2); %Kronecker Delta
    for i=1:2
    for j=1:2
    for k=1:2
    for l=1:2
        C(i,j,k,l) = lambda_o*Finv(j,i)*Finv(l,k) + ...
            mu_o*S(i,k)*S(j,l)-(lambda_o*log(J)...
            - mu_o)*Finv(j,k)*Finv(l,i);
    end
    end
    end
    end   
end