%--------------------------------------------------------------------------
% 2-D Plane Stress Element Routine
% Author: David Vasko
% How to use:
%   Basically just change the element type(elem_typ) and the number of 
%   quadrature points (gauss_pts) and run it to see the output. Curently
%   the output is set to tell you whether or not sepcific testes were
%   passed; give you the defromation gradient at each gauss point; tell you
%   the norm of the error between numerical and analytical values and the
%   rank of teh internal nodal force array and the rank of teh stiffness
%   matrix. It currently plots the deformed and undeformed element.
%   elem_type can be:
%       'lin' for linear element
%       'quad' for quadratic element
%   gauss_pts can be:
%       1 or 3
%   types of deformation: (last argument of gen_deform)
%       'none' - no deformation
%       'scal' - a random pure scaling
%       'rand_corn' - random deformation of just corners
%       'rand_mid' - random deformation for just the middle odes of a
%           quadratic element
%       'rand_all' - random deformation for all nodes in an element. This
%           is the most commonly used type of deformation. Can be used for 
%           both linear adn quadratic elements
%--------------------------------------------------------------------------
%Start Fresh
    close all;clear all;clc;
%RUN 2_D Plane stress element on various deformations
    %first set up our points of our element
        elem_typ = 'quad';
        gauss_pts = 3;
         Xref = [-.5 0; .5 0; 0 1];
        [Xref, xdis] = gen_deform(Xref,elem_typ,'rand_all');
    %Run on this deformation
        planestress_2d(Xref,xdis,elem_typ,gauss_pts);


    
