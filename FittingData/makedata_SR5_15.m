FileName='SR5_15.txt';
Cb = zeros(11,11);%Qb/Qc,AbAc
Cb = Get2DTable(FileName, 6,16);
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cb = Cb.*(AbAc'.^-2*QbQc.^2);
Cb = Cb';

save('SR5_15.mat','Cb','QbQc','AbAc');