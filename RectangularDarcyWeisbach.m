function varargout = RectangularDarcyWeisbach(query, varargin)
if ischar(query)
    query = {query};
end
varargout = cell(1,nargout);
[varargout{:}] = DuctNetwork.RectangularDarcyWeisbachHaaland (query, varargin{:});
end

function varargout = RectangularDarcyWeisbachHaaland(query, varargin)
if ischar(query)
    n=1; query = {query};
elseif iscell(query)
    n = length(query);
else
    varargout = {};return
end
varargout = cell(1,n);
for ii = 1:nargout
    switch query{ii}
        case 'Pdrop'
            q = varargin{1};
            s = varargin{2};
            L = s(1);
            H = s(2);
            W = s(3);
            Dh = 1.3*H^0.625*W^0.625/(H+W)^0.25;
            rho = s(4);
            e = s(5);
            nu = s(6);
          % Re = abs(V)*Dh/nu;
            Re = abs(q)*Dh/H/W/nu;
            
            lambda = 1/(1+exp(-(Re-3750)/250));
            Cf_lam = 64/Re;
            A = (e/3.7/Dh)^3.33;
            B = (6.9/Re)^3;
            T = log10(A+B);
            Cf_turb = (-0.6*T)^(-2);
            Cf = (1-lambda)*Cf_lam + lambda*Cf_turb;
          % M = L/Dh/2*rho*V*abs(V);
            M = L*rho*q*abs(q)/Dh/2/H^2/W^2;
            dP = Cf*M;
            varargout{ii}=dP;
        case 'dPdQ'
            dCf_lamdabsq = -Cf_lam/abs(q);
            dCf_turbdAB = -2*Cf_turb/T/log(10)/(A+B);
            dCf_turbdabsq = -dCf_turbdAB*3*B/abs(q);
            dlambdadabsq = lambda*(1-lambda)/250*Re/abs(q);
            dCfdabsq = (Cf_turb-Cf_lam)*dlambdadabsq + (1-lambda)*dCf_lamdabsq + lambda*dCf_turbdabsq;
            dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
            varargout{ii}=dPdq;
        case 'dPdS'
            dPdL = dP/L;
            Kh = 0.625/H-0.25/(H+W);
          % dDhdH = Dh*Kh;
            dRedH  = Re*Kh;
            dCf_lamdH = -Cf_lam*Kh;
            dAdH = -3.33*A*Kh;
            dBdH = -3*B*Kh;
            dCf_turbdH = dCf_turbdAB*(dAdH+dBdH);
            dlambdadH = lambda*(1-lambda)/250*dRedH;
            dCfdH = (Cf_turb-Cf_lam)*dlambdadH + (1-lambda)*dCf_lamdH + lambda*dCf_turbdH;
            dMdH = M*(0.25/(H+W)-2.625/H);
            dPdH = dCfdH*M+Cf*dMdH;
            Kw = 0.625/W-0.25/(H+W);
          % dDhdW = Dh*Kw;
            dRedW  = Re*Kw;
            dCf_lamdW = -Cf_lam*Kw;
            dAdW = -3.33*A*Kw;
            dBdW = -3*B*Kw;
            dCf_turbdW = dCf_turbdAB*(dAdW+dBdW);
            dlambdadW = lambda*(1-lambda)/250*dRedW;
            dCfdW = (Cf_turb-Cf_lam)*dlambdadW + (1-lambda)*dCf_lamdW + lambda*dCf_turbdW;
            dMdW = M*(0.25/(H+W)-2.625/W);
            dPdW = dCfdW*M+Cf*dMdW;
            dPdrho = dP/rho;
            dCf_turbde = dCf_turbdAB*3.33*A/e;
            dPde = M*lambda*dCf_turbde;
            dCf_lamdnu = Cf_lam/nu;
            dCf_turbdnu = dCf_turbdAB*3*B/nu;
            dlambdadnu = -lambda*(1-lambda)/250*Re/nu;
            dCfdnu = (Cf_turb-Cf_lam)*dlambdadnu + (1-lambda)*dCf_lamdnu + lambda*dCf_turbdnu;
            dPdnu = M*dCfdnu;
            varargout{ii}=[dPdL, dPdH, dPdW, dPdrho, dPde, dPdnu];
        case 'Model_Description'
            varargout{ii}='Rectangular Straight Duct Using Darcy Weisbach Equation by Haaland Approximation';
        case 'Is_Junction'
            varargout{ii}=false;
        case 'Get_Branches'
            varargout{ii}={@DuctNetwork.RectangularDarcyWeisbach};
        case 'Branch_Assignment'
            varargout{ii}={1};
        case 'Parameter_Assignment'
            varargout{ii}={1:6};
        case 'Parameter_Description'
            varargout{ii}={'Length(m)','Height(m)','Width(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
        case 'Is_Shared_Parameter'
            varargout{ii}=[false,false,false,true,true,true];
        case 'Is_Identified_Parameter'
            varargout{ii}=[1,0,0,0,0,0];
        otherwise
            varargout{ii}=[];
    end
end
end

function varargout = RectangularDarcyWeisbachChurchill(query, varargin)
if ischar(query)
    n=1; query = {query};
elseif iscell(query)
    n = length(query);
else
    varargout = {};return
end
varargout = cell(1,n);
for ii = 1:n
    switch query{ii}
        case 'Pdrop'
            q = varargin{1};
            s = varargin{2};
            L = s(1);
            H = s(2);
            W = s(3);
            Dh = 1.3*H^0.625*W^0.625/(H+W)^0.25;
            rho = s(4);
            e = s(5);
            nu = s(6);
          % Re = abs(V)*Dh/nu;
            Re = abs(q)*Dh/H/W/nu;
            
            T1 = power((7/Re),0.9);
            T2 = T1 +(0.27*e/Dh);
            T3 = -2.457*log(T2);
            A = T3^16;
            B = power((37530/Re),16);
            T4 = power((8/Re),12);
            T5 = power(A+B,-1.5);
            Cf = 8*power(T4+T5,1/12);
          % M = L/Dh/2*rho*V*abs(V);
            M = L*rho*q*abs(q)/Dh/2/H^2/W^2;
            dP = Cf*M;
            varargout{ii}=dP;
        case 'dPdQ'
            dAdabsq = 35.3808*A/T3/T2*T1/abs(q);
            dBdabsq = -16*B/abs(q);
            dT4dabsq = -12*T4/abs(q);
            dT5dabsq = -1.5*T5/(A+B)*(dAdabsq+dBdabsq);
            dCfdabsq = 1/12*Cf/(T4+T5)*(dT4dabsq+dT5dabsq);
            dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
            varargout{ii}=dPdq;
        case 'dPdS'
            dPdL = dP/L;
            Kh = 0.625/H-0.25/(H+W);
          % dDhdH = Dh*Kh;
            dRedH  = Re*Kh;
            dT4dH = -12*T4*Kh;
            dT2dH = -0.9*T1*Kh+0.27*e/Dh*(0.25/(H+W)-0.625/H);
            dAdH = -2.457*16*A/T3/T2*dT2dH;
            dBdH = -16*B*Kh;
            dT5dH = -1.5*T5/(A+B)*(dAdH+dBdH);
            dCfdH = 1/12*Cf/(T4+T5)*(dT4dH+dT5dH);
            dMdH = M*(0.25/(H+W)-2.625/H);
            dPdH = dCfdH*M+Cf*dMdH;
            Kw = 0.625/H-0.25/(H+W);
          % dDhdW = Dh*Kw;
            dRedW  = Re*Kw;
            dT4dW = -12*T4*Kw;
            dT2dW = -0.9*T1*Kw+0.27*e/Dh*(0.25/(H+W)-0.625/W);
            dAdW = -2.457*16*A/T3/T2*dT2dW;
            dBdW = -16*B*Kw;
            dT5dW = -1.5*T5/(A+B)*(dAdW+dBdW);
            dCfdW = 1/12*Cf/(T4+T5)*(dT4dW+dT5dW);
            dMdW = M*(0.25/(H+W)-2.625/W);
            dPdW = dCfdW*M+Cf*dMdW;
            dPdrho = dP/rho;
            % dAde = -2.457*16*0.27*A/T3/T2/D;
            % dT5de = -1.5*T5/(A+B)*dAde;
            % dCfde = 1/12*Cf/(T4+T5)*dT5de;
            % dPde = M*dCfde;
            dPde = 1.32678*dP*T3^15/T2/(T4+T5)/(A+B)^2.5/D;
            dT4dnu = 12*T4/nu;
            dAdnu = -2.457*16*0.9*A/T3/T2*T1/nu;
            dBdnu = 16*B/nu;
            dT5dnu = -1.5*T5/(A+B)*(dAdnu+dBdnu);
            dCfdnu = 1/12*Cf/(T4+T5)*(dT4dnu+dT5dnu);
            dPdnu = M*dCfdnu;
            varargout{ii}=[dPdL, dPdH, dPdW, dPdrho, dPde, dPdnu];
        case 'Model_Description'
            varargout{ii}='Rectangular Straight Duct Using Darcy Weisbach Equation by Churchill Approximation';
        case 'Is_Junction'
            varargout{ii}=false;
        case 'Get_Branches'
            varargout{ii}={@DuctNetwork.RectangularDarcyWeisbach};
        case 'Branch_Assignment'
            varargout{ii}={1};
        case 'Parameter_Assignment'
            varargout{ii}={1:6};
        case 'Parameter_Description'
            varargout{ii}={'Length(m)','Height(m)','Width(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
        case 'Is_Shared_Parameter'
            varargout{ii}=[false,false,false,true,true,true];
        case 'Is_Identified_Parameter'
            varargout{ii}=[1,0,0,0,0,0];
        otherwise
            varargout{ii}=[];
    end
end
end
