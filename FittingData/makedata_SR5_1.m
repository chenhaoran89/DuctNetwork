Cb = zeros(11,3,3);%Qb/Qc,AbAc,AsAc
Cs = zeros(11,3,3);%Qs/Qc,AbAc,AsAc

FileName='SR5_1.txt';
for ii=1:3
    Cb(:,:,ii) = Get3DTable(FileName, (ii-1)*3+6,(ii)*3+5)';
    Cs(:,:,ii) = Get3DTable(FileName, (ii-1)*3+21,(ii)*3+20)';
end

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.25      0.5      1.0];
AsAc = [0.5      0.75      1.0];

[N1,N2,N3] = size(Cb);
Cb = Cb.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,N2,1),N1,1,N3).^2;

[N1,N2,N3] = size(Cs);
Cs = Cs.*repmat(reshape(QsQc',N1,1,1),1,N2,N3)./repmat(reshape(AsAc,1,1,N3),N1,N2,1).^2;

save('SR5_1.mat','Cb','Cs','QbQc','QsQc','AsAc','AbAc');