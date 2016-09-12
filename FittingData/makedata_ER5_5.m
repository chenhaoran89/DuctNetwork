FileName='ER5_5.txt';
Cb = zeros(4,11);%Qb/Qc,AbAc
Cb = Get2DTable(FileName, 6,9);
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.25,0.5,0.75,1];
Cb = Cb.*(AbAc'.^-2*QbQc.^2);
Cb = Cb';

save('ER5_5.mat','Cb','QbQc','AbAc');