diretorio = 'C:\physiobank\database\edb\extracted\'; 
list = dir(diretorio);
V1 = [];
cont =0;
records = java_array('java.lang.String',11);
for index = 3:size(list)
    fileName = list(index).name;
  
   if(not(isequal('e0403.mat', fileName)))
        strcat(diretorio,fileName);
        matObj = matfile(strcat(diretorio,fileName));
        varlist = who(matObj);

        for index2 = 1:size(varlist)
            varName = varlist{index2};
            if(isequal('V1',varName))
                V1 = [V1, matObj.(varName)];
                cont = cont +1;
                records(cont) = java.lang.String(fileName);
            end
        end 
    end
    
end
cont

V1Rocha = zeros(1,16);
V1RochaNormal = zeros(1,16);
V1RochaIsch = zeros(1,16);
[size1, size2] = size(V1);
for index = 1:size2
    V1Rocha = cat(1, V1Rocha, V1(index).Datasets.Rocha);
end
[size3, size4] = size( V1Rocha);
for index2 = 1:size3
    x = V1Rocha(index2);
    t = V1Rocha(index2,15);
    if(t == 0)
        V1RochaNormal = cat(2, V1RochaNormal, x);
    else
        V1RochaIsch = cat(2, V1RochaIsch, x);
    end
end

[sizeN1, sizeN2] = size(V1RochaNormal);
[sizeI1, sizeI2] = size(V1RochaIsch);

sizeN2+sizeI2
save v1RochaSet.mat V1Rocha;
save v1RochaNormalSet.mat V1RochaNormal;
save v1RochaIschSet.mat V1RochaIsch;
    