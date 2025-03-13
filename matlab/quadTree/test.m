%%

dx = [1 ; 2 ; 1];
dy = 2*ones(4,1);
dz = [2 ; 4; 6];
Nza = 0;

coarseGrid = Grid3D(dx,dy,dz,Nza);
Kmax = 1;

obj = quadtreeGrid(coarseGrid,Kmax);

Klayer = [0 1 0]';
obj.setCellsLayered(Klayer);

%%
obj.setNodes

obj.setEdges

obj.setFaces

obj.setCellVolume

obj.setEdgeLength
    %%
    temp = zeros(size(obj.xEdges));
    i1 = 1; i2 = obj.nEdges(1);
    temp(obj.indE(i1:i2)) = obj.E(i1:i2)
    %%
    temp = zeros(size(obj.yEdges));
    i1 = i2+1; i2 = i2+obj.nEdges(2);
    temp(obj.indE(i1:i2)) = obj.E(i1:i2)
    %%
    temp = zeros(size(obj.zEdges));
    i1 = i2+1; i2 = i2+obj.nEdges(3);
    temp(obj.indE(i1:i2)) = obj.E(i1:i2)
    
%%
obj.setFaceArea
    %%
    temp = zeros(size(obj.xFaces));
    i1 = 1; i2 = obj.nFaces(1);
    temp(obj.indF(i1:i2)) = obj.F(i1:i2)
    %%
    temp = zeros(size(obj.yFaces));
    i1 = i2+1; i2 = i2+obj.nFaces(2);
    temp(obj.indF(i1:i2)) = obj.F(i1:i2)
    %%
    temp = zeros(size(obj.zFaces));
    i1 = i2+1; i2 = i2+obj.nFaces(3);
    temp(obj.indF(i1:i2)) = obj.F(i1:i2)
    
%%
obj.setDiv
%%
obj.setCurl
%%
obj.setDualEdgeLength
%%
obj.setDualFaceArea