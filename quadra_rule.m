function [PandW] = quadra_rule(pts)
%Pand w will be of form [theta1 theta2 weight] for each row of PandW
%dont need theta3 becaue theta3 = 1 - theta1 - tehta2
    switch(pts)
        case 1
            PandW = [1/3 1/3 1];
        case 3
            PandW = [1/2 1/2 1/3; 0 1/2 1/3; 1/2 0 1/3];
        otherwise
            disp('ERROR: Gauss points must be 2 or 3');
            PandW = [];
            return;
    end
end