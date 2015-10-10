%Declare Constants that will be used in neohookean eqns
    lambda_o = 3;
    mu_o = 2;
%Right Cauchy Green Deformation Tensor
%and other usefull constants
    RCGDT = F.'*F;
    J = det(F);
    I_1 = trace(RCGDT);
    Finv = F^-1;
    lnJ = log(J);
%strain energy density   
    CiJkL2 = zeros(3,3,3,3);
    CIJKL2 = CiJkL2;
    S = eye(3); %Kronecker Delta
    for i=1:3
    for j=1:3
    for k=1:3
    for l=1:3
       % CiJkL2(i,j,k,l) = lambda_o*Finv(j,i)*Finv(l,k) + mu_o*S(i,k)*S(j,l);
    %Original tangent moduli
        CiJkL2(i,j,k,l) = lambda_o*Finv(j,i)*Finv(l,k) + ...
            mu_o*S(i,k)*S(j,l)-(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i);
    %Get the altered form ofteh tangent moduli (based on second piola)
        CIJKL2(i,j,k,l) = 0.5*Finv(j,i)*Finv(l,k)*(CiJkL2(i,j,k,l)-S(i,k)*S(j,l));
                         
    end
    end
    end
    end 
CiJkL = zeros(3,3,3,3);
CIJKL = CiJkL;
for k=1:3
for l=1:3
CiJkL(:,:,k,l) = lambda_o*Finv'*Finv(l,k) + mu_o*S(:,k)*S(:,l)'...
    -(lambda_o*lnJ - mu_o)*Finv(l,:)'*Finv(:,k)';
CIJKL(:,:,k,l) = 0.5*Finv'*Finv(l,k).*(CiJkL(:,:,k,l)-S(:,k)*S(:,l)');
end
end
% CiJkL(:,:,2,1) = lambda_o*Finv'*Finv(1,2) + mu_o*S(:,2)*S(:,1)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,3,1) = lambda_o*Finv'*Finv(1,3) + mu_o*S(:,3)*S(:,1)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,1,2) = lambda_o*Finv'*Finv(2,1) + mu_o*S(:,1)*S(:,2)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,2,2) = lambda_o*Finv'*Finv(2,2) + mu_o*S(:,2)*S(:,2)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,3,2) = lambda_o*Finv'*Finv(2,3) + mu_o*S(:,3)*S(:,2)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,1,3) = lambda_o*Finv'*Finv(3,1) + mu_o*S(:,1)*S(:,3)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,2,3) = lambda_o*Finv'*Finv(3,2) + mu_o*S(:,2)*S(:,3)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
% CiJkL(:,:,3,3) = lambda_o*Finv'*Finv(3,3) + mu_o*S(:,3)*S(:,3)'...
%     -(lambda_o*lnJ - mu_o)*Finv(j,k)*Finv(l,i)';
%=-------------------------------------------------------------------------
[~,gi] = get_duals(G_i,g_i);
        [~,~,CIJKL] = neoHookean(G_i,g_i);
        Cijkl = zeros(3,3,3,3);
        for i=1:3
        for I=1:3
        for j=1:3
        for J=1:3
        for k=1:3
        for K=1:3
        for l=1:3
        for L=1:3
            Cijkl(i,j,k,l) = Cijkl(i,j,k,l) +...
                CIJKL(I,J,K,L)*gi(I,i)*gi(J,j)*gi(K,k)*gi(L,l);
        end
        end
        end
        end
        end
        end
        end
        end
        
        Cijkl2 = zeros(3,3,3,3);

        for k=1:3
        for l=1:3
        for J=1:3
        for I=1:3
            tx = reshape(CIJKL(I,J,:,:),3,3).*(gi(:,k)*gi(:,l)');
            T = sum(sum(tx));
            M = (gi(I,:)'*gi(J,:));               
            Cijkl2(:,:,k,l) = Cijkl2(:,:,k,l) + M*T;
        end
        end
        end
        end
