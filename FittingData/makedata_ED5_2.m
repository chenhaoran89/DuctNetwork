Cb = zeros(11,11,11);%Qb/Qc,AbAc,AsAc
Cs = zeros(11,11,11);%Qs/Qc,AbAc,AsAc

FileName='ED5_2.txt';
for ii=1:11
    Cb(:,:,ii) = Get3DTable(FileName, (ii-1)*11+6,(ii)*11+5)';
    Cs(:,:,ii) = Get3DTable(FileName, (ii-1)*11+133,(ii)*11+132)';
end

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];

[N1,N2,N3] = size(Cb);
Cb = Cb.*repmat(reshape(QbQc',N1,1,1),1,N2,N3)./repmat(reshape(AbAc,1,N2,1),N1,1,N3).^2;

[N1,N2,N3] = size(Cs);
Cs = Cs.*repmat(reshape(QsQc',N1,1,1),1,N2,N3)./repmat(reshape(AsAc,1,1,N3),N1,N2,1).^2;

save('ED5_2.mat','Cb','Cs','QbQc','QsQc','AsAc','AbAc');