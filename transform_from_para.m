function [x,y,z] = transform_from_para(X,para_pt,elem_typ)
%decidewhat shape function to use
    switch elem_typ
        case 'lin'
            shapefunc = @linTri_Na;
            num_nodes = 3;
        case 'quad'
            shapefunc = @quadTri_Na; 
            num_nodes = 6;
        otherwise
            disp('ERROR: we are using linear(lin) or quadratic(quad)');
            x = 0; y = 0;
            return;
    end
%Nodal X and Y values   
    xia = X(1,:);
    yia = X(2,:);
    zia = X(3,:);
%Gaussian x and y values
        x = 0; y = 0; z = 0;
        Na = shapefunc(para_pt(1),para_pt(2));
        for a=1:num_nodes
            x = x + xia(a)*Na(a);
            y = y + yia(a)*Na(a);
            z = z + zia(a)*Na(a);
        end
end