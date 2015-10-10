function [Xref, xdis] = gen_deform(Xref,elem_typ,deform_typ,amt)
switch deform_typ
    case 'scal'
        xdis = Xref*(1+rand()*amt);
        [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,false);
    case 'none'
        xdis = Xref;
        [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,false);
    case 'rand_corn'
        xdis = Xref + 2*(rand(3,3)-0.5*ones(3,3))*amt;
        [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,false);
    case 'rand_mid'
        xdis = Xref;
        [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,true);
    case 'rand_all'
        xdis = Xref + 2*(rand(3,3)-0.5*ones(3,3))*amt;
        [Xref,xdis] = quad_midpts(Xref,xdis,elem_typ,true);
    otherwise
        disp('invalid deform type');
end
end
