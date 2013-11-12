lista = dir('files.c')
for i=1:length(lista)    
    lista(i).name
    mex (lista(i).name);
end
