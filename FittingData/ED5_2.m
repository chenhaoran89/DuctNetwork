function [dP, dPdQ, dPdS]=ED5_2(q, s, Selection)
% q = [Qs,Qb,Qc];
% s = [Ds,Db,Dc,rho];
ModifiedTable = 0;
if (ModifiedTable == 1)
    gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
    dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
    dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
    switch Selection
        case 'b'
            GridVec = {DuctNetwork.Table_ED5_2.QbQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
            ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
            dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
            dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            GridVec = {DuctNetwork.Table_ED5_2.QsQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
            ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
            dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
            dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
    end
else
    switch Selection
        case 'b'
            gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
            dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
            dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
            GridVec = {DuctNetwork.Table_ED5_2.QbQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
            ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
            dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
            dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
            dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
            dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
            GridVec = {DuctNetwork.Table_ED5_2.QsQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
            ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
            dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
            dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
    end
end
end