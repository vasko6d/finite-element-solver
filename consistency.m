function [errors] = consistency(Xref, xdis,gauss_pts,elem_typ,TRI,h)
%declare teh preterbation h
    dispher = false;
    if nargin <=5;
        h = 10^-5;
        dispher = true;
    end
%decide if doing for a single element or a mesh
    switch length(TRI(:,1))
        case 1
            only_one = true;
        otherwise
            only_one = false;
    end
if only_one
%Get values for analytical values and delcare limits
    Pext = [0;0;0];
    [~,f,~,K] = WfK(Xref,xdis,gauss_pts,elem_typ,Pext);
    num_nodes = length(Xref(1,:));
    ind_max=3;
%preterb each node by a small bit
        fh = zeros(ind_max,num_nodes);
            for a=1:num_nodes
                for i=1:ind_max            
                %reset x plus h and x minus h
                    xph = xdis;
                    xmh = xdis;
                %preterb the node in teh X direcction
                    xph(i,a) = xph(i,a) + h;
                    xmh(i,a) = xmh(i,a) - h;
                    Wph = WfK(Xref,xph,gauss_pts,elem_typ,Pext,true);
                    Wmh = WfK(Xref,xmh,gauss_pts,elem_typ,Pext,true);
                    fh(i,a) = (Wph - Wmh)/(2*h);
                end
            end
%             disp('numerical f'); disp(fh); disp('analytical f'); disp(f);
            ef = max(max(abs(f-fh)));
            if dispher , disp('error between f and fh');disp(ef); end
            %e_s(z) = norm(f-fh,'fro');
%preterb each node by a small bit
        Kh = zeros(ind_max,num_nodes,ind_max,num_nodes);
            for a=1:num_nodes
            for b=1:num_nodes
                for i=1:ind_max
                for k=1:ind_max
                %reset x plus h and x minus h
                    xph = xdis;
                    xmh = xdis;
                %preterb the node in teh X direcction
                    xph(k,b) = xph(k,b) + h;
                    xmh(k,b) = xmh(k,b) - h;
                    [~,fph,~,~] = WfK(Xref,xph,gauss_pts,elem_typ,Pext,true);
                    [~,fmh,~,~] = WfK(Xref,xmh,gauss_pts,elem_typ,Pext,true);
                    Kh(i,a,k,b) = (fph(i,a) - fmh(i,a))/(2*h);            
                end
                end
            end
            end            
%             disp('numerical K'); disp(unwrap_K(Kh));
%             disp('analytical K'); disp(unwrap_K(K));
            eK = max(max(max(max(abs(K-Kh)))));
            if dispher, disp('error between K and Kh');disp(eK); end
%output the errors to user
    errors = [ef eK];
else
    disp('hello');
%Get values for analytical values and delcare limits
    Pext = [0;0;0];
    [~,F,~,K] = assemble(Xref,xdis,TRI,Pext,gauss_pts,elem_typ);
    num_nodes = length(Xref(1,:));
    ind_max=3;
%preterb each node by a small bit
        fh = zeros(ind_max,num_nodes);
            for a=1:num_nodes
                for i=1:ind_max            
                %reset x plus h and x minus h
                    xph = xdis;
                    xmh = xdis;
                %preterb the node in teh X direcction
                    xph(i,a) = xph(i,a) + h;
                    xmh(i,a) = xmh(i,a) - h;
                    Wph = assemble(Xref,xph,TRI,Pext,...
                        gauss_pts,elem_typ);
                    Wmh = assemble(Xref,xmh,TRI,Pext,...
                        gauss_pts,elem_typ);
                    fh(i,a) = (Wph - Wmh)/(2*h);
                end
            end
%             disp('numerical f'); disp(fh); disp('analytical f'); disp(f);
            ef = max(max(abs(F-fh)));
            if dispher , disp('error between f and fh');disp(ef); end
            %e_s(z) = norm(f-fh,'fro');
%preterb each node by a small bit
        Kh = zeros(ind_max,num_nodes,ind_max,num_nodes);
            for a=1:num_nodes
            for b=1:num_nodes
                for i=1:ind_max
                for k=1:ind_max
                %reset x plus h and x minus h
                    xph = xdis;
                    xmh = xdis;
                %preterb the node in teh X direcction
                    xph(k,b) = xph(k,b) + h;
                    xmh(k,b) = xmh(k,b) - h;
                    [~,fph,~,~] = assemble(Xref,xph,TRI,Pext,...
                        gauss_pts,elem_typ);
                    [~,fmh,~,~] = assemble(Xref,xmh,TRI,Pext,...
                        gauss_pts,elem_typ);
                    Kh(i,a,k,b) = (fph(i,a) - fmh(i,a))/(2*h);            
                end
                end
            end
            end            
%             disp('numerical K'); disp(unwrap_K(Kh));
%             disp('analytical K'); disp(unwrap_K(K));
            eK = max(max(max(max(abs(K-Kh)))));
            if dispher, disp('error between K and Kh');disp(eK); end
%output the errors to user
    errors = [ef eK];
end
end