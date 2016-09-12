function [dP, dPdQ, dPdS]=ER5_3(q, s, Selection)
% q=[qs,qb,qc];
% s=[Hs,Ws,Hb,Wb,rho];
gExp =  0.5*s(5)*(q(3)/(s(1)*s(2)))^2;
dgdq =  [0,0,s(5)*q(3)/(s(1)*s(2))^2];
dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
switch Selection
    case 'b'
        Cb_Table = DuctNetwork.Table_ER5_3.Cb;
        GridVec = {DuctNetwork.Table_ER5_3.QbQc};
        ZExp = [q(2)/q(3)];
        dZdq = [0,1/q(3),-q(2)/q(3)^2];
        dZds = [0,0,0,0,0];
        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cb_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
    case 's'
        Cs_Table = DuctNetwork.Table_ER5_3.Cs;
        GridVec = {DuctNetwork.Table_ER5_3.QsQc};
        ZExp = [q(1)/q(3)];
        dZdq = [1/q(3),0,-q(1)/q(3)^2];
        dZds = [0,0,0,0,0];
        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cs_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
    otherwise
        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
end
end

