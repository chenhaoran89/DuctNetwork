function [f,dfdq,dfds] = Interp_Gradient(GridVector, InterpTable, ZExp, dZdq, dZds, gExp, dgdq, dgds, q,s)
% Z = ZExp(q,s)
% C = Interpolation(GridVector,InterpTable,Z)
% g = gExp(q,s)
% f = g*C;
% dfdq = g*dCdZ*dZdq + C*dgdq;
% dfds = g*dCdZ*dZds + C*dgds;
Z = reshape(ZExp(q,s),1,[]);
NZ = length(Z);
StepSize = cell(1,NZ);IndexInTable = cell(1,NZ);ZMesh = cell(1,NZ);IndexInMesh = cell(1,NZ);
for ii = 1:NZ
    if Z(ii)>GridVector{ii}(end)
        IndexInTable{ii}=length(GridVector{ii})*[1;1];
        StepSize{ii} = 2*(Z(ii)-GridVector{ii}(end));
        ZMesh{ii}=[GridVector{ii}(end);Z(ii);2*Z(ii)-GridVector{ii}(end)];
    elseif Z(ii)<GridVector{ii}(1)
        IndexInTable{ii}=[1;1];
        StepSize{ii} = -2*(Z(ii)-GridVector{ii}(1));
        ZMesh{ii}=[2*Z(ii)-GridVector{ii}(1);Z(ii);GridVector{ii}(1)];
    else
        IndexInTable{ii}=find(Z(ii)<GridVector{ii},1)*[1;1]+[-1;0];
        StepSize{ii} = range(GridVector{ii}(IndexInTable{ii}));
        ZMesh{ii}=[GridVector{ii}(IndexInTable{ii}(1));GridVector{ii}(IndexInTable{ii}(end))];
    end
    IndexInMesh{ii}=[1;2];
end
CfMesh = InterpTable(IndexInTable{:});
lambda = cellfun(@(v,x)(x-v(1))/(v(2)-v(1)),ZMesh,num2cell(Z));
for ii=1:NZ
    Index_1 = IndexInMesh; Index_1{ii}=1;
    Index_2 = IndexInMesh; Index_2{ii}=1;
    Index_3 = IndexInMesh; Index_3{ii}=1;
    CfMesh(Index_3{:}) = CfMesh(Index_2{:});
    CfMesh(Index_2{:}) = CfMesh(Index_1{:})*lambda(ii)+CfMesh(Index_3{:})*(1-lambda(ii));
    IndexInMesh{ii}=[1;2;3];
end
Index_mid = num2cell(2*ones(1,NZ));
Cf=CfMesh(Index_mid{ii});
cell_Grad = cell(1,NZ);
[cell_Grad{:}] = gradient(CfMesh,StepSize{:});
dCfdZ = cellfun(@(M)M(Index_mid{:}),cell_Grad);
g = gExp(q,s);
dZdq = dZdq(q,s);
dZds = dZds(q,s);
dgdq = dgdq(q,s);
dgds = dgds(q,s);
f = Cf*g;
dfdq = g*dCfdZ*dZdq + Cf*dgdq;
dfds = g*dCfdZ*dZds + Cf*dgds;
end