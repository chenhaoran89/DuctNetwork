function [dP, dPdQ, dPdS]=ER5_5(q, s, Selection)
gExp =  0.5*s(5)*(q(3)/(s(1)*s(4)))^2;
dgdq =  [0,0,s(5)*q(3)/(s(1)*s(4))^2];
dgds =  [-2/s(1),0,0,-2/s(4),1/s(5)]*gExp;
switch Selection
    case 'b'
        GridVec = {DuctNetwork.Table_ER5_5.QbQc,DuctNetwork.Table_ER5_5.AbAc};
        ZExp = [q(2)/q(3),s(2)/s(4)];
        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
        dZds = [0,0,0,0,0;0,1/s(4),0,-s(2)/s(4)^2,0];
        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ER5_5.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
    otherwise
        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
end
end
