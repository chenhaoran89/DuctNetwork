duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)','Roughness(mm)','Dynamic Viscosity(m^2/s)'},[1.204,0.15,15.11e-6],[2,3,4]);
duct.AddBranch(0,1);
duct.AddFitting(1,'FanQuadratic',[20,0.2]);
duct.AddBranch(1,0);
duct.AddBranch(1,0);
duct.AddFitting([1,2,3],'SR5_13_Junction',rand(1,6));
% duct.S = [1.20400000000000;0.150000000000000;1.51100000000000e-05;20;0.456424600851747;0.713795583233135;0.884405045275995;0.720855670816932;0.0186127747263861;0.674776467128286];
[X,Q,P]=duct.Sim();
Q
% %%
% h=1e-8;
% %X0 = [0.88;0.22];
% X0 = rand(2,1);
% 
% Branch = 1;
% [dP,dPdQ,dPdS]=duct.BranchPressureDrop(Branch,X0,duct.S);
% jac2 = DuctNetwork.Jacobian(@(x)duct.BranchPressureDrop(Branch,x,duct.S),X0,1e-6);
% jac3 = DuctNetwork.Jacobian(@(x)duct.BranchPressureDrop(Branch,x,duct.S),X0,-1e-6);
% [dPdQ;jac2;jac3]
