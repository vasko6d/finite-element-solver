function [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,rand_mid)
if strcmp(elem_typ,'quad')
    Xtemp = zeros(6,3);
    Xtemp(1:3,:) = Xref;
    Xtemp(4,:) = mean(Xtemp(1:2,:));
    Xtemp(5,:) = mean(Xtemp(2:3,:));
    Xtemp(6,:) = mean(Xtemp(1:2:3,:));
    xtemp = zeros(6,3);
    xtemp(1:3,:) = xdis;
    xtemp(4,:)=mean(xtemp(1:2,:));
    xtemp(5,:)=mean(xtemp(2:3,:));
    xtemp(6,:)=mean(xtemp(1:2:3,:));
    Xref = Xtemp;
    xdis = xtemp;
    if rand_mid
        xdis(4,:)=xdis(4,:)+2*((rand(1,3)-0.5*ones(1,3)))*.05;
        xdis(5,:)=xdis(5,:)+2*((rand(1,3)-0.5*ones(1,3)))*.05;
        xdis(6,:)=xdis(6,:)+2*((rand(1,3)-0.5*ones(1,3)))*.05;
    end
end
end