function [G_i,g_i]= get_tans(X,x,Na,Na_a)
%set up repeted shape function matrixes to make sums easy
    I3 = eye(3);
    Nr_mat = zeros(3,length(Na)*3);
    Ns_mat = zeros(3,length(Na)*3);
    for i=1:length(Na)
        ind = 1+(i-1)*3;
        Nr_mat(:,ind:ind+2) = I3*Na_a(i,1);
        Ns_mat(:,ind:ind+2) = I3*Na_a(i,2);
    end  
%Put the the nodal positions in verticle form
    n = length(X(1,:));
    X_line = reshape(X,n*3,1);
    x_line = reshape(x,n*3,1);
%Get the reference tangent vectors
    Gr = Nr_mat*X_line;
    Gs = Ns_mat*X_line;    
    G3 = cross(Gr,Gs);    
%get the deformed tangent vectors
    gr = Nr_mat*x_line;
    gs = Ns_mat*x_line;    
    g3 = cross(gr,gs); 
%Get the magnitude of gr X gs and return everything
    sqrtA = [norm(G3),norm(g3)];
    G_i = [Gr,Gs,G3/sqrtA(1)];
    g_i= [gr,gs,g3/sqrtA(2)];
end