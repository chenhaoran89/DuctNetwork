duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)'},1.204);
duct.AddNode('Node')
duct.AddBranch(0,1);% Branch 1
duct.AddBranch(1,0);% Branch 2
duct.AddBranch(1,0);% Branch 3
duct.AddFitting(1,'PressureSource',00)
duct.AddFitting(-2,'PressureSource',50)
duct.AddFitting(-3,'PressureSource',150)
duct.AddFitting([1,2,3],'RectangularYJunction',[0.4,0.5,0.3,0.3])
duct.Sim();
