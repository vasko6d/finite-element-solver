function [Na, Na_a] = linTri_Na(r,s)
%Get the Shape Functions
    N1 = r;
    N2 = s;
    N3 = 1 -r -s;
    Na = [N1 N2 N3];
%Compute The Derivatives of Na
    N1_1 = 1;
    N1_2 = 0;
    N2_1 = 0;
    N2_2 = 1;
    N3_1 = -1;
    N3_2 = -1;
    Na_a = [N1_1 N1_2; N2_1 N2_2; N3_1 N3_2];
end