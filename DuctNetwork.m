classdef DuctNetwork < handle
    % DuctNetwork implements all algorithms for balancing duct network
    
    properties
        n; % number of nodes without ground
        n_NodeDescription; % cell array of size n by 1, text description on each node
        P; % array of size n by 1, pressure vector on each node, Pa
        b; % number of branches in the network
        b_BranchDescription; % cell array of size b by 1, text description on each branch
        Q; % array of size b by 1, flow rate vector on each branch, m^3/s
        A; % Associate matrix of size n by b
        % A(i,j)==1 means branch j leaves node i
        % A(i,j)==0 means branch j is not connected with node i
        % A(i,j)==-1 means branch j enters node i
        t; % dimension of null space
        U; % null space basis of matrix A, of size b by t, AU=0, U'U=1
        b_Pdrop; %cell array of size b by 1 for the pressure drop in terms of q and s
        %if there exist multiple Pdrop functions, use a cell array to store each of them
        %only record the original Pdrop function for the fitting, handle the fitting direction in pressure drop calculation
        b_dPdQ; %cell array of size b by 1 for the partial pressure drop over q in terms of q and s,
        %only record the original dPdQ function for the fitting, handle the fitting direction in pressure drop calculation
        b_dPdS; %cell array of size b by 1 for the partial pressure drop over s in terms of q and s
        %only record the original dPdS function for the fitting, handle the fitting direction in pressure drop calculation
        b_Qidx; %cell array of size b by 1 for the index of dependent Q of Pdrop, dPdQ and dPdS functions, so use Q(abs(b_Qidx{b})).*sign(b_Qidx{b}) in these functions.
        % Qidx{ii}[jj] can be negative, the sign represents the fitting direction w.r.t. branch direction
        b_Sidx; %cell array of size b by 1 for the index of dependent S of Pdrop, dPdQ and dPdS functions, so use S(b_Sidx{b}) in these functions.
        s; % number of parameters in the model
        s_ParamDescription; % cell array of size b by 1, text description on each parameter
        S; % array of size s by 1, parameter vector for whole system
        s_m; % array of size s by 1, s_m(i)=0 means no need to identify S(i), s_m(i)>0 is the multiplicity of this parameter that need to be identified
        s_MultiS; % cell array of size s by 1, containing values of multiplicities of each parameter in row vectors, Use [s_MultiS{:}]' to obtain all in one column vector
        
    end
    
    methods
        function obj = DuctNetwork()
            obj.n = 0;
            obj.n_NodeDescription = cell(0,1);
            obj.P = zeros(0,1);
            obj.b = 0;
            obj.b_BranchDescription = cell(0,1);
            obj.Q = zeros(0,1);
            obj.A = zeros(0,0);
            obj.t = 0;
            obj.U = zeros(0,0);
            obj.b_Pdrop = cell(0,1);
            obj.b_dPdQ = cell(0,1);
            obj.b_dPdS = cell(0,1);
            obj.b_Qidx = cell(0,1);
            obj.b_Sidx = cell(0,1);
            obj.s = 0;
            obj.s_ParamDescription = cell(0,1);
            obj.S = zeros(0,1);
            obj.s_m = zeros(0,1);
            obj.s_MultiS = cell(0,1);
        end
        
        branch = ParamDependency(para_idx); % this functions querry the dependent branches for given parameter index.

        function AddBranch(obj,FromNode,ToNode,varargin)
            [~, b_Curr]=size(obj.A);
            Branch_Idx = b_Curr+1;
            function SetNode(Node,a)
                if ischar(Node)
                    switch lower(Node)
                        case {'atm','gnd'}
                            Node=0;
                        otherwise
                            Node = find(cellfun(@(S)strcmp(S,Node), obj.n_NodeDescription),1);
                            if isempty(Node)
                                Node=0;
                            end
                    end
                end
                if Node>0
                    obj.A(Node,Branch_Idx) = a;
                end
            end
            SetNode(FromNode,1);
            SetNode(ToNode,-1);
            [obj.n,obj.b]=size(obj.A);
            obj.U = null(obj.A);
            obj.t = size(obj.U,2);
            obj.b_Pdrop{Branch_Idx,1}=cell(0,1);
            obj.b_dPdQ{Branch_Idx,1}=cell(0,1);
            obj.b_dPdS{Branch_Idx,1}=cell(0,1);
            obj.b_Qidx{Branch_Idx,1}=cell(0,1);
            obj.b_Sidx{Branch_Idx,1}=cell(0,1);
            obj.b_BranchDescription{Branch_Idx,1}=cell(0,1);
            if nargin>3
                obj.AddFitting(Branch_Idx,varargin{:});
            end
        end
        
        function Node=AddNode(obj,varargin)
            Node = obj.n+1;
            obj.n = Node;
            if nargin>1
                obj.n_NodeDescription{Node,1} = varargin{1};
            end
        end
        
        function AddFitting(obj,Branch,Model,Param) %Branch has direction, positve means fitting direction is same as the branch direction, negative means opposite
            if isa(Model,'function_handle')
                if Model('Is_Junction')
                    % separate the junction into three/four junction branches, and add each of them separately, parameters are shared by all of them.
                    CenterNode=find(arrayfun(@(i)isempty(setdiff(find(obj.A(i,:)),abs(Branch))),1:obj.n),1);
                    if isempty(CenterNode)
                        disp('no such junction exist, fail to create junction'); return
                    end    
                    Branch = -sign(obj.A(CenterNode,abs([1,2,3]))).*abs(Branch);
                    branch_function_handle = Model('Get_Branches');
                    branch_assignment = Model('Branch_Assignment');
                    param_assignment = Model('Parameter_Assignment');
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(Param));
                    for ii = 1:length(ParamDescription)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(cellfun(@(S)strcmp(S,ParamDescription{ii}),obj.s_ParamDescription),1);
                            if isempty(s_Idx) % shared parameter is first included in the model
                                obj.s = obj.s+1; %add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; %add the parameter at the end
                            obj.s_ParamDescription{obj.s} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    for ii = 1:length(branch_function_handle)
                        BranchList = Branch(branch_assignment{ii});
                        PrimaryBranch = abs(BranchList(1));
                        obj.b_BranchDescription{PrimaryBranch} = [obj.b_BranchDescription{PrimaryBranch};{feval(branch_function_handle{ii},'Model_Description')}];
                        obj.b_Pdrop{PrimaryBranch} = [obj.b_Pdrop{PrimaryBranch};{@(q,s)feval(branch_function_handle{ii},{'Pdrop','dPdQ','dPdS'},q,s)}];
                        obj.b_dPdQ{PrimaryBranch} = [obj.b_dPdQ{PrimaryBranch};{@(q,s)feval(branch_function_handle{ii},'dPdQ',q,s)}];
                        obj.b_dPdS{PrimaryBranch} = [obj.b_dPdS{PrimaryBranch};{@(q,s)feval(branch_function_handle{ii},'dPdS',q,s)}];
                        obj.b_Qidx{PrimaryBranch} = [obj.b_Qidx{PrimaryBranch};{BranchList}];
                        obj.b_Sidx{PrimaryBranch} = [obj.b_Sidx{PrimaryBranch};{Param_Idx(param_assignment{ii})}];
                    end
                    obj.S(Param_Idx) = Param;
                    obj.s_m(Param_Idx) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param,Param_Idx,'UniformOutput',false);
                else
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(Param));
                    for ii = 1:length(bShared)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(cellfun(@(S)strcmp(S,ParamDescription{ii}),obj.s_ParamDescription),1);
                            if isempty(s_Idx)
                                obj.s = obj.s+1; %add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; %add the parameter at the end
                            obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    
                    BranchList = Branch;
                    PrimaryBranch = abs(Branch(1));
                    obj.b_BranchDescription{PrimaryBranch} = [obj.b_BranchDescription{PrimaryBranch};{Model('Model_Description')}];
                    obj.b_Pdrop{PrimaryBranch} = [obj.b_Pdrop{PrimaryBranch};{@(q,s)Model({'Pdrop','dPdQ','dPdS'},q,s)}];
                    obj.b_dPdQ{PrimaryBranch} = [obj.b_dPdQ{PrimaryBranch};{@(q,s)Model('dPdQ',q,s)}];
                    obj.b_dPdS{PrimaryBranch} = [obj.b_dPdS{PrimaryBranch};{@(q,s)Model('dPdS',q,s)}];
                    obj.b_Qidx{PrimaryBranch} = [obj.b_Qidx{PrimaryBranch};{BranchList}];
                    obj.b_Sidx{PrimaryBranch} = [obj.b_Sidx{PrimaryBranch};{Param_Idx}];
                    
                    obj.S(Param_Idx) = Param;
                    obj.s_m(Param_Idx) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param,obj.b_Sidx{Branch}{end},'UniformOutput',false);
                    
                end
            elseif ischar(Model)
                Model = DuctNetwork.FittingDatabase(Model);
                obj.AddFitting(Branch,Model,Param);
            elseif iscell(Model)
                cellfun(@(a,b)obj.AddFitting(Branch,a,b),Model,Param);
            end
        end

        function [dP, dPdX, dPdS] = BranchPressureDrop(obj,Branch_Idx,Q,S)
            N = length(obj.b_BranchDescription{Branch_Idx});
            dP = zeros(N,1);dPdQ = zeros(N,obj.b);dPdS = zeros(N,obj.s);dPdX = zeros(N,obj.t);
            for k = 1:N
                dir = sign(obj.b_Qidx{Branch_Idx}{k});
                idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                idx_S = obj.b_Sidx{Branch_Idx}{k};
                [dP(k), dPdQ(k,idx_Q), dPdS(k,idx_S)] = feval(obj.b_Pdrop{Branch_Idx}{k},Q(idx_Q),S(idx_S));
                dP(k) = dir(1)*dP(k);
                dPdX(k,:) = (dir(1)*dPdQ(k,idx_Q).*dir)*obj.U(idx_Q,:);
                dPdS(k,:) = dir(1)*dPdS(k,:);
            end
            dP = sum(dP,1);
            dPdX = sum(dPdX,1);
            dPdS = sum(dPdS,1);
        end
            
        function [e, dedX, dedS] = res_StateEquation(obj,X,S)
            % Input:
            % X is obj.t dimensional vector of internal state
            % S is obj.s dimensional vector of system parameter used in simulation
            % Output:
            % e is the residual of State Equations of size obj.t by 1
            % dedX is the Jacobian of e wrt internal state X of size obj.t by obj.t
            % dedS is the Jacobian of e wrt system parameter S of size obj.t by obj.s
            X = real(X);
            [dP, dPdX, dPdS]=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S),(1:obj.b)','UniformOutput',false);
            dP = cell2mat(dP); dPdX = cell2mat(dPdX); dPdS = cell2mat(dPdS);
            e = obj.U'*dP;
            dedX = obj.U'*dPdX;
            dedS = obj.U'*dPdS;
        end
        
        function [X,Q,P] = Sim(obj,Param_Value, Param_Description)
            S_value = obj.S;
            if nargin>=3
                S_value(cellfun(@(Str) find(strcmp(Str,obj.s_ParamDescription)),Param_Description)) = Param_Value;
            elseif nargin>=2
                S_value = Param_Value;
            end
            
            X0 = rand(obj.t,1);
            options = optimoptions(@fsolve,'Display','none',...
                'Algorithm','trust-region-reflective',...
                'TolFun',1e-12,'TolX',1e-12,...
                'MaxIterations',obj.t*5,...
                'Jacobian','on','DerivativeCheck','off');
            [X,~,exitflag] = fsolve(@(x) obj.res_StateEquation(x,S_value),X0,options);
            if exitflag<=0
                disp('Change to numerical Jacobian!')
                options = optimoptions(options,'Jacobian','off');
                X = fsolve(@(x) obj.res_StateEquation(x,S_value),X0,options);
            end
            dP=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S_value),(1:obj.b)','UniformOutput',false);
            dP=cell2mat(dP);
            Q = obj.U*X;
            P = (obj.A*obj.A')\obj.A*dP;
            obj.Q = Q;
            obj.P = P;
        end
        
    end
    
    methods (Static = true)
        function ModelFcn = FittingDatabase(ModelStr)
            ModelFcn=str2func(['DuctNetwork.',ModelStr]);
        end
        
        function varargout = FanQuadratic(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Fan Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.FanQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Max Pressure(Pa)','Max Flow(m^3/s)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,1];
                    case 'Pdrop'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        if q<0
                            varargout{ii} = s(1);
                        else
                            varargout{ii} = s(1)*(q^2/s(2)^2-1);
                        end
                    case 'dPdQ'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        if q<0
                            varargout{ii} = 0;
                        else
                            varargout{ii} = s(1)*2*q/s(2)^2;
                        end
                    case 'dPdS'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        if q<0
                            varargout{ii} = [1,0];
                        else
                            varargout{ii} = [q^2/s(2)^2-1, -s(1)*2*q^2/s(2)^3];
                        end
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = DuctQuadratic (query, varargin)
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
                    case 'Model_Description'
                        varargout{ii}='Duct Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.DuctQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'Resistance(Pa/(m^3/s)^2)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    case 'Pdrop'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        varargout{ii} = s(1)*q*abs(q);
                    case 'dPdQ'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        varargout{ii} = s(1)*2*abs(q);
                    case 'dPdS'
                        if (~exist('q','var')), q = varargin{1};s = varargin{2}; end;
                        varargout{ii} = q*abs(q);
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularDarcyWeisbach( query, varargin)
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
                    case 'Model_Description'
                        varargout{ii}='Circular Straight Duct Using Darcy Weisbach Equation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Dynamic Viscosity(m^2/s)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0];
                    case 'Pdrop'
                        if (~exist('dP','var')), DataInitiation(); end;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        if (~exist('dP','var')), DataInitiation(); end;
                        dRedq = D/nu/Area;
                        invRedRedq=dRedq/Re;
                        dAdq = 35.3808*T1^15/T2*T21*invRedRedq;
                        dBdq = -16*B*invRedRedq;
                        dCfdq = Cf/(T3+T4)*(-T3*invRedRedq-T4/8/(A+B)*(dAdq+dBdq));
                        %dPdq = Cf*L/D*rho*V/Area + L/D*rho/2*V^2*dCfdq;
                        dPdq = T5*(Cf/Area+V/2*dCfdq);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        if (~exist('dP','var')), DataInitiation(); end;
                        S1 = -2.457*16*T1^15/T2;
                        dAdD = S1*(0.9*T21-0.27*e/D)/D;
                        dBdD = 16*B/D;
                        
                        G = Cf/12/(T3+T4);
                        H1 = 12*T3;
                        H2 = 1.5/(A+B)^2.5;
                        
                        dPdrho = dP/rho;
                        dPdL = dP/L;
                        
                        dCfdD = G*(H1/D-H2*(dAdD+dBdD));
                        dPdD = dP*(-3/D+dCfdD/Cf);
                        dAdnu = S1*(0.9*T21/nu);
                        dBdnu = 16*B/nu;
                        dCfdnu = G*(H1/nu-H2*(dAdnu+dBdnu));
                        dPdnu = dP/Cf*dCfdnu;
                        
                        dAde = S1*0.27/D;
                        dCfde = -G*H2*dAde;
                        dPde = dP/Cf*dCfde;
                        
                        dPds = [dPdL, dPdD, dPdrho, dPde, dPdnu];
                        varargout{ii}=dPds;
                    otherwise
                        varargout{ii}=[];
                end
            end
            
            function DataInitiation()
                q = varargin{1};
                s = varargin{2};
                L = s(1);
                D = s(2);
                rho = s(3);
                e = s(4);
                nu = s(5);
                
                Area = pi*(D/2)^2;
                V = q/Area;
                Re =abs(V)*D/nu;
                
                
                T21 = realpow((7/Re),0.9);
                T2 = T21 +(0.27*e/D);
                T1 = -2.457*reallog(T2);
                A = T1^16;
                
                B = realpow((37530/Re),16);
                T3 = realpow((8/Re),12);
                T4 = 1/realpow(A+B,1.5);
                Cf = 8*realpow(T3+T4,1/12);
                
                T5 = L/D*rho*abs(V);
                %dP = Cf*L/D*rho/2*V^2;
                dP = Cf*T5/2*V;
                if ~isreal(dP)
                    disp(dP)
                end
            end
        end
        
        function varargout = CircularTJunction( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Main,@DuctNetwork.CircularTJunction_Main,@DuctNetwork.CircularTJunction_Branch};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4];[2,1,3,4];[3,1,2,4]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(4,1);
                end
            end
        end
        
        function varargout = CircularTJunction_Main( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Main of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Main};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    case 'Pdrop'
                        if (~exist('dP','var')), DataInitiation(); end;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        if (~exist('dPdQ','var')), DataInitiation(); end;
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        if (~exist('dPdS','var')), DataInitiation(); end;
                        varargout{ii}=dPdS;
                end
            end
            function DataInitiation()
                q = reshape(varargin{1},1,[]);
                s = reshape(varargin{2},1,[]);
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SD5-18, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.Calc_SD5_18(abs(q([1,2,3])), s([1,2,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.Calc_SD5_9(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ED5-3 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SD5-9 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ED5-3 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.Calc_ED5_3(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ED5-4, use Cb
                        if s(1)>=s(2) % D1>D2, use Cb1
                            [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.Calc_ED5_4(abs(q([1,2,3])), s([1,2,3,4]), 'b1');
                        else % D1<D2, use Cb2
                            [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_ED5_4(abs(q([2,1,3])), s([2,1,3,4]), 'b2');
                        end
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible                        
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularTJunction_Branch( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Branch of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Branch};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Diameter(m)','Main 1 Diameter(m)','Main 2 Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    case 'Pdrop'
                        if (~exist('dP','var')), DataInitiation(); end;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        if (~exist('dPdQ','var')), DataInitiation(); end;
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        if (~exist('dPdS','var')), DataInitiation(); end;
                        varargout{ii}=dPdS;
                end
            end
            function DataInitiation()
                q = reshape(varargin{1},1,[]);
                s = reshape(varargin{2},1,[]);
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_SD5_9(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.Calc_SD5_9(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ED5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SD5-18 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.Calc_ED5_3(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the main 1, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_ED5_3(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible                        
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
                
            end
        end
    end
    
    properties (Constant)
        ED5_3 = load('FittingData/ED5_3.mat');
        SD5_9 = load('FittingData/SD5_9.mat');
        ED5_4 = load('FittingData/ED5_4.mat');
        SD5_18 = load('FittingData/SD5_18.mat');
    end
    
    methods (Static = true)
        function [dP, dPdQ, dPdS]=Calc_ED5_3(q,s,Selection)
%             syms q1 q2 q3 s1 s2 s3 s4 real;
%             VarList = [q1,q2,q3,s1,s2,s3,s4];
%             A = pi*s3^2/4;V = q3/A; gExp = 0.5*s4*V^2;
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,2/q(3)]*gExp(q,s);
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b'
                    if s(3)<=0.25 %Dc<= 0.25m
                        Cb_Table = DuctNetwork.ED5_3.Cb_part1;
                    else
                        Cb_Table = DuctNetwork.ED5_3.Cb_part2;
                    end
                    Interpolant = griddedInterpolant({DuctNetwork.ED5_3.QbQc,DuctNetwork.ED5_3.AbAc,DuctNetwork.ED5_3.AsAc},Cb_Table,'linear','nearest');
%                     ZExp = [q2/q3;(s2/s3)^2;(s1/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(InterpCb, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                case 's'
                    if s(3)<=0.25 %Dc<= 0.25m
                        Cs_Table = DuctNetwork.ED5_3.Cs_part1;
                    else
                        Cs_Table = DuctNetwork.ED5_3.Cs_part2;
                    end
                    Interpolant = griddedInterpolant({DuctNetwork.ED5_3.QsQc,DuctNetwork.ED5_3.AbAc,DuctNetwork.ED5_3.AsAc},Cs_Table,'linear','nearest');
%                     ZExp = [q1/q3;(s2/s3)^2;(s1/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(InterpCs, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_SD5_9(q,s,Selection)
%             syms q1 q2 q3 s1 s2 s3 s4 real;
%             VarList = [q1,q2,q3,s1,s2,s3,s4];
%             A = pi*s3^2/4;V = q3/A; gExp = 0.5*s4*V^2;
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,2/q(3)]*gExp(q,s);
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b'
                    Interpolant = griddedInterpolant({DuctNetwork.SD5_9.QbQc,DuctNetwork.SD5_9.AbAc},DuctNetwork.SD5_9.Cb,'linear','nearest');
%                     ZExp = [q2/q3;(s2/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(2)/q(3);(s(2)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                case 's'
                    Interpolant = griddedInterpolant({DuctNetwork.SD5_9.QsQc,DuctNetwork.SD5_9.AsAc},DuctNetwork.SD5_9.Cs,'linear','nearest');
%                     ZExp = [q1/q3;(s1/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_ED5_4(q,s,Selection)
%             syms q1 q2 q3 s1 s2 s3 s4 real;
%             VarList = [q1,q2,q3,s1,s2,s3,s4];
%             A = pi*s3^2/4;V = q3/A; gExp = 0.5*s4*V^2;
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,2/q(3)]*gExp(q,s);
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b1'
                    Interpolant = griddedInterpolant({DuctNetwork.ED5_4.QbQc,DuctNetwork.ED5_4.AbAc,DuctNetwork.ED5_4.AbAc},DuctNetwork.ED5_4.Cb1,'linear','nearest');
%                     ZExp = [q1/q3;(s1/s3)^2;(s2/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2;(s(2)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                case 'b2'
                    Interpolant = griddedInterpolant({DuctNetwork.ED5_4.QbQc,DuctNetwork.ED5_4.AbAc,DuctNetwork.ED5_4.AbAc},DuctNetwork.ED5_4.Cb2,'linear','nearest');
%                     ZExp = [q2/q3;(s1/s3)^2;(s2/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(2)/q(3);(s(1)/s(3))^2;(s(2)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_SD5_18(q,s,Selection)
%             syms q1 q2 q3 s1 s2 s3 s4 real;
%             VarList = [q1,q2,q3,s1,s2,s3,s4];
%             A = pi*s3^2/4;V = q3/A; gExp = 0.5*s4*V^2;
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,2/q(3)]*gExp(q,s);
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b'
                    Interpolant = griddedInterpolant({DuctNetwork.SD5_18.QbQc,DuctNetwork.SD5_18.AbAc},DuctNetwork.SD5_18.Cb,'linear','nearest');
%                     ZExp = [q1/q3;(s1/s3)^2];
%                     [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, gExp, VarList, q, s, 0.01);
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.NumericalEval_InterpolantGradient(Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s, 0.01);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [f,dfdq,dfds] = NumericalEval_InterpolantGradient(Coeff_Interpolant, ZExp, dZdq, dZds, gExp, dgdq, dgds, q,s, h)
            % X = [q,s]
            % Z = ZExp(q,s)
            % C = Coeff_Interpolant(Z)
            % g = gExp(q,s)
            % f = C*g
            % dfdq = g*dCdZ*dZdq + C*dgdq;
            % dfds = g*dCdZ*dZds + C*dgds;
            
            Z0 = reshape(feval(ZExp,q,s),1,[]);
            NZ = length(Z0);
            Z = num2cell([Z0+h/2;Z0;Z0-h/2],1);
            Y = Coeff_Interpolant(Z);
            GRAD = cell(1,NZ);
            [GRAD{:}] = gradient(Y,h/2);
            Idx = num2cell(2*ones(1,NZ));
            C = Y(Idx{:});
            dCdZ = cellfun(@(M)M(Idx{:}),GRAD,'UniformOutput',true);
            g = feval(gExp,q,s);
            dZdq = feval(dZdq,q,s);
            dZds = feval(dZds,q,s);
            dgdq = feval(dgdq,q,s);
            dgds = feval(dgds,q,s);
            f = C*g;
            dfdq = g*dCdZ*dZdq + C*dgdq;
            dfds = g*dCdZ*dZds + C*dgds;
            
%             X0 = [reshape(q,1,[]),reshape(s,1,[])];
%             Nq = length(q);
%             Z0=double(subs(ZExp,VarList,X0))';
%             NZ = length(Z0);
%             Z = num2cell([Z0+h/2;Z0;Z0-h/2],1);
%             Y = Coeff_Interpolant(Z);
%             GRAD = cell(1,NZ);
%             [GRAD{:}] = gradient(Y,h/2);
%             Idx = num2cell(2*ones(1,NZ));
%             C = Y(Idx{:});
%             dCdZ = cellfun(@(M)M(Idx{:}),GRAD,'UniformOutput',true);
%             g = double(subs(gExp,VarList,X0));
%             dZdX = double(subs(jacobian(ZExp),VarList,X0));
%             dgdX = double(subs(jacobian(gExp),VarList,X0));
%             f = C*g;
%             T_ZExp = cell2mat(arrayfun(@(m)has(symvar(ZExp)',m),VarList,'UniformOutput',false));
%             T_gExp = cell2mat(arrayfun(@(m)has(symvar(gExp)',m),VarList,'UniformOutput',false));
%             dfdX = g*dCdZ*dZdX*T_ZExp + C*dgdX*T_gExp;
%             dfdq = dfdX(:,1:Nq);
%             dfds = dfdX(:,Nq+1:end);
        end
    end
end

