classdef DuctNetwork < handle
    % DuctNetwork implements all algorithms for balancing duct network 
    
    properties
        n; % number of nodes without ground
        n_NodeDescription; % cell array of size n by 1, text description on each node
        P; % array of size n by 1, pressure vector on each node, Pa
        b; % number of branches in the network
        b_BranchDescription; % cell array of size b by 1, text description on each branch
        Q; % array of size b by 1, flow rate vector on each branch, m^3/s
        nb_A; % Associate matrix of size n by b 
            % A(i,j)==1 means branch j leaves node i
            % A(i,j)==0 means branch j is not connected with node i
            % A(i,j)==-1 means branch j enters node i
        b_Pdrop; %cell array of size b by 1 for the pressure drop in terms of q and s
        %if there exist multiple Pdrop functions, use a cell array to store each of them
        b_dPdQ; %cell array of size b by 1 for the partial pressure drop over q in terms of q and s, 
        b_dPdS; %cell array of size b by 1 for the partial pressure drop over s in terms of q and s
        b_Qidx; %cell array of size b by 1 for the index of dependent Q of Pdrop, dPdQ and dPdS functions, so use Q(b_Qidx{b}) in these functions.
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
            obj.n_NodeDescription = {};
            obj.P = [];
            obj.b = 0;
            obj.b_BranchDescription = {};
            obj.Q = [];
            obj.nb_A = [];
            obj.b_Pdrop = {};
            obj.b_dPdQ = {};
            obj.b_dPdS = {};
            obj.b_Qidx = {};
            obj.b_Sidx = {};
            obj.s = 0;
            obj.s_ParamDescription = {};
            obj.S = [];
            obj.s_m = [];
            obj.s_MultiS = {};
            
        end
        
        branch = ParamDependency(para_idx); % this functions querry the dependent branches for given parameter index.
    end
    
end

