function [Na, Na_a] = quadTri_Na(r,s)
%Write out the lambdas
    lambda1 = r;
    lambda2 = s;
    lambda3 = 1 -r -s;
%Get the Shape Functions
    N4 = 4*lambda1*lambda2;
    N5 = 4*lambda2*lambda3;
    N6 = 4*lambda1*lambda3;
    N1 = lambda1*(2*lambda1-1);
    N2 = lambda2*(2*lambda2-1);
    N3 = lambda3*(2*lambda3-1);
    Na = [N1 N2 N3 N4 N5 N6];
%Compute The Derivatives of Na wrt r and s (lambda1 and lambda2)
    N1_r = 4*lambda1-1;
    N1_s = 0;
    N2_r = 0;
    N2_s = 4*lambda2-1;
    N3_r = 4*(lambda1 + lambda2)-3;
    N3_s = 4*(lambda1 + lambda2)-3;
    N4_r = 4*lambda2;
    N4_s = 4*lambda1;
    N5_r = -4*lambda2;
    N5_s = 4-4*lambda1-8*lambda2;
    N6_r = 4-4*lambda2-8*lambda1;
    N6_s = -4*lambda1;
    Na_a = [N1_r N1_s; N2_r N2_s; N3_r N3_s; 
            N4_r N4_s; N5_r N5_s; N6_r N6_s];
end