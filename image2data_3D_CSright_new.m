% rootdir = 'C:\Users\44797\Desktop\bedforms4\';
% caseNum = 3;
% simulationNum = 4;
% simname = sprintf('%sRomain_%d\\Simulation_%d',rootdir, caseNum, simulationNum);
% foldername = sprintf('%sRomain_%d\\Simulation_%d\\cs_y_base',rootdir, caseNum, simulationNum);

csstep = 1;
zaspect = 0.25;

simname = outputFolder_simulation;
foldername = outputFolder_simulation;

cd(foldername);   


% H = 235;
% ImageDis = 995;
% INum = 200;

tic;
colourMatrixR = nan(H,ImageDis,INum);
colourMatrixG = nan(H,ImageDis,INum);
colourMatrixB = nan(H,ImageDis,INum);
sliceTemp = 0;

colourMatrixFacies = nan(H,ImageDis,INum);


for j = 1 : INum

    cd(foldername);
    
    fileName = ['y_', num2str(j),'.png'];

    I = imread(fileName); 

    R = I(:,:,1);
    G = I(:,:,2);
    B = I(:,:,3);

    sliceTemp = sliceTemp + 1;
    colourMatrixR(:,:,sliceTemp) = R;
    colourMatrixG(:,:,sliceTemp) = G;
    colourMatrixB(:,:,sliceTemp) = B;
    
    
    F = nan(H,ImageDis);
    
%     [r1,c1] = find(R==165); 
%     for i = 1 : numel(r1)
%         F(r1(i),c1(i)) = 1;
%     end
%     
%     [r2,c2] = find(R==246);
%     for i = 1 : numel(r2)
%         F(r2(i),c2(i)) = 2;
%     end
%     
%     [r3,c3] = find(G==190);
%     for i = 1 : numel(r3)
%         F(r3(i),c3(i)) = 3;
%     end
%     
%     [r4,c4] = find(G==236);
%     for i = 1 : numel(r4)
%         F(r4(i),c4(i)) = 4;
%     end
%     
%     [r6,c6] = find(R==234);
%     for i = 1 : numel(r6)
%         F(r6(i),c6(i)) = 5;
%     end
%     
%     [r5,c5] = find(R==255);
%     for i = 1 : numel(r5)
%         F(r5(i),c5(i)) = NaN;
%     end



    % R G B
    % 96 64 44
    % 139 115 114
    % 254 190 0
    % 246 170 110
    % 254 236 0




    [r1,c1] = find(R==96); 
    for i = 1 : numel(r1)
        F(r1(i),c1(i)) = 1; 
    end
    
    [r2,c2] = find(R==139);
    for i = 1 : numel(r2)
        F(r2(i),c2(i)) = 2;
    end
    
    [r3,c3] = find(G==190);
    for i = 1 : numel(r3)
        F(r3(i),c3(i)) = 3;
    end
    
    [r4,c4] = find(G==170);
    for i = 1 : numel(r4)
        F(r4(i),c4(i)) = 4;
    end
    
    [r6,c6] = find(G==236);
    for i = 1 : numel(r6)
        F(r6(i),c6(i)) = 5;
    end

    [r7,c7] = find(G==111);
    for i = 1 : numel(r7)
        F(r7(i),c7(i)) = 6;
    end
    
    [r5,c5] = find(R==255);
    for i = 1 : numel(r5)
        F(r5(i),c5(i)) = NaN;
    end

    colourMatrixFacies(:,:,sliceTemp) = F;
    
        

end


toc;
disp('New matrix completed');

%%

CS_YZ = permute(colourMatrixFacies,[1 3 2]); % Y-by-Z array
CS_XY = permute(colourMatrixFacies,[3 2 1]); % X-by-Y array
CS_XZ = permute(colourMatrixFacies,[2 3 1]); % X-by-Z array

cd(simname);
load('TopoINfO.mat');
load('ContourLineINfO.mat');
ZBED_yz = permute(ZBEDDraw,[1 3 2]); % Y-by-Z array
ZBED_xy = permute(ZBEDDraw,[3 2 1]); % X-by-Y array
ZBED_xz = permute(ZBEDDraw,[2 3 1]); % X-by-Z array


GridSize_y = 100;
resolution_y = GridSize_y/size(CS_YZ,2);
CenterShift_y = round(GridSize_y/2);
axisThis_y = -CenterShift_y + resolution_y : resolution_y : GridSize_y - CenterShift_y;

GridSize_x = 100;
resolution_x = GridSize_x/size(CS_YZ,3);
CenterShift_x = round(GridSize_x/2);
axisThis_x = -CenterShift_x + resolution_x : resolution_x : GridSize_x - CenterShift_x;


figymin = min(min(min(ZBED_yz)));
figymax = max(max(max(ZBED_yz)));
resolution_z = (figymax-figymin)/size(CS_YZ,1);
CenterShift_z = round((figymax-figymin)/2);
axisThis_z = figymin+resolution_z : resolution_z : figymax;

[xGrid,yGrid] = meshgrid(axisThis_y, axisThis_z);
[xGridhori,yGridhori] = meshgrid(axisThis_x, axisThis_y);
[xGrid_x,yGrid_z] = meshgrid(axisThis_x, axisThis_z);


%% cross section along x
% SurfGridSpace = 0.5; % 1,0.5
% trimEdge = 10;
trimNum = trimEdge/SurfGridSpace;

figxmin = min(min(XSurf)) + trimEdge;
figxmax = max(max(XSurf)) - trimEdge;
figymin = min(min(YSurf)) + trimEdge;
figymax = max(max(YSurf)) - trimEdge;
figzmin = min(min(min(ZBEDDraw)));
% figzmin = -12;
figzmax = max(max(max(ZBEDDraw)));

ZBEDDrawTop = ZBEDDraw(:,:,1);
linesW = 0.25;

zthisH = axisThis_z(1,size(CS_XY, 3));
dis = abs(zrefList - zthisH);
[~, I] = min(dis);
thisHorContLines = ContourLinesList{I};

% thiscs = CS_XY(:,:,k);

z3D = zeros(size(xGridhori,1),size(xGridhori,2));
z3D = z3D + zthisH;
    
mymap = [
  96/256, 64/256, 44/256
  139/256, 115/256, 114/256
  254/256, 190/256, 0/256
  246/256, 170/256, 110/256
  254/256, 236/256, 0/256
  37/256, 111/256, 229/256];


fg = figure(1);
  % =============== cross section along x ==================
for j = 1+trimNum : csstep: size(ZBED_xz,3)-trimNum
    zthis = ZBED_xz(:,:,j);
    ythis = YSurf(j,1);

    dis = abs(axisThis_y - ythis);
    [v, I] = min(dis);
    thiscs = colourMatrixFacies(:,:,I);
    thiscs = rot90(thiscs,2);
    thiscs = fliplr(thiscs);


    y3D = zeros(size(xGrid_x,1),size(xGrid_x,2));
    y3D = y3D + ythis;

    zshow = yGrid_z;
    zshow(zshow>zthisH) = NaN;



    set(0,'CurrentFigure',fg);
    set(gcf,'color','w');

    % ================ Surface ================
    ax2 = axes;
    ZBEDDrawTopshow = ZBEDDrawTop;
    ZBEDDrawTopshow(ZBEDDrawTopshow > zthisH + 0.2) = NaN;

    ytemp = YSurf(:,1);
    [rowdel, ~] = find(ytemp < ythis);
    rowNumberStart = max(rowdel+1);
    XSurfnew = XSurf;
    YSurfnew = YSurf;
    ZBEDDrawTopshownew = ZBEDDrawTopshow;

    XSurfnew(rowdel,:) = [];
    YSurfnew(rowdel,:) = [];
    ZBEDDrawTopshownew(rowdel,:) = [];

    Bedformplot = surfl(ax2, XSurfnew, YSurfnew, ZBEDDrawTopshownew,[60,100],[.1,.4,.3,100]);
    axis([figxmin figxmax figymin figymax figzmin figzmax]);
    
    daspect([1 1 zaspect]);


    colormap(ax2,'gray');
    shading interp
    material dull;
    axis off
    hold on

    ZBEDDrawTopLeft = ZBEDDrawTopshow(j:end-trimNum,1+trimNum);
    ZBEDDrawTopLeft(ZBEDDrawTopLeft > zthisH) = zthisH;
    xLeft = XSurf(j: end-trimNum, 1+trimNum);
    yLeft = YSurf(j: end-trimNum, 1+trimNum);

    xbox1=[figxmin figxmin xLeft' figxmin];
    ybox1 = [figymax ythis yLeft' figymax];
    zbox1=[figzmin figzmin ZBEDDrawTopLeft' figzmin];

    [~,cdel] = find(isnan(zbox1));
    zbox1(:,cdel) = [];
    xbox1(:,cdel) = [];
    ybox1(:,cdel) = [];
    fill3(ax2, xbox1,ybox1,zbox1,[.7 .7 .7], 'EdgeColor','k')



    ZBEDDrawTopRight = ZBEDDrawTop(j, 1+trimNum:end-trimNum);
    ZBEDDrawTopRight(ZBEDDrawTopRight > zthisH) = zthisH;
    xRight = XSurf(j, 1+trimNum:end-trimNum);
    yRight = YSurf(j, 1+trimNum:end-trimNum);

    xbox2=[figxmax figxmin xRight figxmax];
    ybox2=[ythis ythis yRight ythis];
    zbox2=[figzmin figzmin ZBEDDrawTopRight figzmin];

    [~,cdel] = find(isnan(zbox2));
    zbox2(:,cdel) = [];
    xbox2(:,cdel) = [];
    ybox2(:,cdel) = [];

    fill3(xbox2,ybox2,zbox2,[.9 .9 .9],'EdgeColor','k')


    xLeftBack = XSurf(1+trimNum:end-trimNum,end-trimNum);
    yLeftBack = YSurf(1+trimNum:end-trimNum,end-trimNum);
    ZBEDDrawTopLeftBack = ZBEDDrawTop(1+trimNum:end-trimNum,end-trimNum);

    xbox3=[figxmax figxmax xLeftBack' figxmax];
    ybox3=[figymax figymin yLeftBack' figymax];
    zbox3=[figzmin figzmin ZBEDDrawTopLeftBack' figzmin];

    [~,cdel] = find(isnan(zbox3));
    zbox3(:,cdel) = [];
    xbox3(:,cdel) = [];
    ybox3(:,cdel) = [];

    fill3(xbox3,ybox3,zbox3,[.7 .7 .7],'EdgeColor','k')

    xbox4=[figxmin figxmax figxmax figxmin];
    ybox4=[figymin figymin figymax figymax];
    zbox4=[figzmin figzmin figzmin figzmin];
    fill3(xbox4,ybox4,zbox4,[.7 .7 .7],'EdgeColor','k')


    %  ================================

  ax1 = axes;
  s = surf(ax1, xGrid_x, y3D, zshow, thiscs, 'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour
  % view([0 90]);
  axis([figxmin figxmax figymin figymax figzmin figzmax]);
  hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
  ax1.Visible = 'off';
  ax1.XTick = [];
  ax1.YTick = [];

  colormap(ax1, mymap);
  caxis manual
  caxis([1 6]);

  daspect([1 1 zaspect]);

  hold on;
  ythis = YSurf(j,:);

  zshow2 = zthis;
  zshow2(zshow2>zthisH) = NaN;

  XSurfnew = XSurf(j,:);
  plot3(ax1, XSurfnew, ythis, zshow2, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);

  % =============== cross section along y ===========
  zthis = ZBED_yz(:,:,1+trimNum);
  xthis = XSurf(1,1+trimNum);
  xGridNew = xGrid*(-1);

  dis = abs(axisThis_x - xthis);
  [~, I] = min(dis);
  thiscs = CS_YZ(:,:,I);
  thiscs = rot90(thiscs,2);

  x3D = zeros(size(yGrid,1),size(yGrid,2));
  x3D = x3D + xthis;

  zshow = yGrid;
  zshow(zshow>zthisH) = NaN;
    
  ytemp = xGridNew(1,:);
  [~, coldel] = find(ytemp < ythis);

  x3D_new = x3D;
  xGridNew_new = xGridNew;
  zshow_new = zshow;
  thiscs_new = thiscs;

  x3D_new(:,coldel) = [];
  xGridNew_new(:,coldel) = [];
  zshow_new(:,coldel) = [];
  thiscs_new(:,coldel) = [];

  surf(ax1, x3D_new, xGridNew_new,zshow_new,thiscs_new,'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour
  % view([0 90]);
  axis([figxmin figxmax figymin figymax figzmin figzmax]);
  colormap(ax1, mymap);
  caxis manual
  caxis([1 6]);


  hold on;
  xshow = XSurf(:,1+trimNum);
  yshow = YSurf(:,1+trimNum);
  zshow2 = zthis;
  zshow2(zshow2>zthisH) = NaN;
  
  [rdel,~] = find(yshow < ythis);
  xshow_new = xshow;
  yshow_new = yshow;
  zshow2_new = zshow2;
  xshow_new(rdel,:) = [];
  yshow_new(rdel,:) = [];
  zshow2_new(rdel,:) = [];

  plot3(ax1, xshow_new, yshow_new, zshow2_new, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);

  view(-38,20);

  hold off;




  % set(gcf,'GraphicsSmoothing','off');
  figname = sprintf('%s\\yy_%s.png', foldername, num2str(j-trimNum));
  print('-dpng','-r600',figname);
  clf;

%     figname2 = sprintf('%s\\z_%s.eps', foldername, num2str(ii-trimNum));
%     print(gcf,'-depsc','-painters',figname2); 
%     epsclean(figname2);   
end
      
      


