load SLAB/slab5.bln
load SLAB/slab10.bln
load SLAB/slab20.bln
load SLAB/slab30.bln
load SLAB/slab40.bln
load SLAB/slab50.bln
load SLAB/slab60.bln
load SLAB/slab70.bln
load SLAB/slab80.bln
slab5  = [slab5; slab5(end,1) 41];
slab10  = [slab10; slab10(end,1) 41];
slab20  = [slab20; slab20(end,1) 41];
slab30  = [slab30; slab30(end,1) 41];
slab40  = [slab40; slab40(end,1) 41];
slab50  = [slab50; slab50(end,1) 41];
slab60  = [slab60; slab60(end,1) 41];
slab70  = [slab70; slab70(end,1) 41];
slab80  = [slab80; slab80(end,1) 41];
slab100 = slab80;
slab100(:,1) = slab100(:,1)+.295*2;
%slab120 = slab100;
%slab120(:,1) = slab120(:,1)+.295*2;
%slab140 = slab120;
%slab140(:,1) = slab140(:,1)+.295*2;
%slab160 = slab140;
%slab160(:,1) = slab160(:,1)+.295*2;
slabs{1} = slab5;
slabs{2} = slab10;
slabs{3} = slab20;
slabs{4} = slab30;
slabs{5} = slab40;
slabs{6} = slab50;
slabs{7} = slab60;
slabs{8} = slab70;
slabs{9} = slab80;
slabs{10} = slab100;
%slabs{11} = slab120;
%slabs{12} = slab140;
%slabs{13} = slab160;
depths = [5,10,20,30,40,50,60,70,80, 100]% 120 140 160];
