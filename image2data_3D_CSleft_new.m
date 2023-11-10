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
    
    %fileName = ['y_', num2str(1),'.png'];
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


%% cross section along y
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
linesW = 0.1;

zthisH = axisThis_z(1,size(CS_XY, 3));
dis = abs(zrefList - zthisH);
[~, I] = min(dis);
thisHorContLines = ContourLinesList{I};

% thiscs = CS_XY(:,:,k);

z3D = zeros(size(xGridhori,1),size(xGridhori,2));
z3D = z3D + zthisH;

% R G B
% 96 64 44
% 139 115 114
% 254 190 0
% 246 170 110
% 254 236 0
% 37 111 229

mymap = [
  96/256, 64/256, 44/256
  139/256, 115/256, 114/256
  254/256, 190/256, 0/256
  246/256, 170/256, 110/256
  254/256, 236/256, 0/256
  37/256, 111/256, 229/256];
 

 
 fg = figure(1);
 for ii = 1+trimNum : csstep : size(ZBED_yz,3)-trimNum
 % for ii = 1+trimNum : 1+trimNum
   
    zthis = ZBED_yz(:,:,ii);
    xthis = XSurf(1,ii);
    xGridNew = xGrid*(-1);

    dis = abs(axisThis_x - xthis);
    [~, I] = min(dis);
    thiscs = CS_YZ(:,:,I);
    thiscs = rot90(thiscs,2);

    x3D = zeros(size(yGrid,1),size(yGrid,2));
    x3D = x3D + xthis;

    zshow = yGrid;
    zshow(zshow>zthisH) = NaN;
      
      
    set(0,'CurrentFigure',fg);
    set(gcf,'color','w');

    % ================ Surface ================


    ax2 = axes;
    ZBEDDrawTopshow = ZBEDDrawTop;
    ZBEDDrawTopshow(ZBEDDrawTopshow > zthisH + 0.2) = NaN;
    
    xtemp = XSurf(1,:);
    [~, coldel] = find(xtemp < xthis);
    colNumberStart = max(coldel+1);
    XSurfnew = XSurf;
    YSurfnew = YSurf;
    ZBEDDrawTopshownew = ZBEDDrawTopshow;
    
    XSurfnew(:,coldel) = [];
    YSurfnew(:,coldel) = [];
    ZBEDDrawTopshownew(:,coldel) = [];

    Bedformplot = surfl(ax2, XSurfnew, YSurfnew, ZBEDDrawTopshownew,[60,100],[.1,.4,.3,100]);
    axis([figxmin figxmax figymin figymax figzmin figzmax]);


    daspect([1 1 zaspect]);

    colormap(ax2,'gray');
    shading interp
    material dull;
    axis off
    hold on

    ZBEDDrawTopLeft = ZBEDDrawTopshownew(1+trimNum:end-trimNum,1);
    ZBEDDrawTopLeft(ZBEDDrawTopLeft > zthisH) = zthisH;

    
    xLeft = XSurf(1+trimNum: end-trimNum, ii);
    yLeft = YSurf(1+trimNum: end-trimNum, ii);

    xbox1=[xthis xthis xLeft' xthis];
    ybox1 = [figymax figymin yLeft' figymax];
    zbox1=[figzmin figzmin ZBEDDrawTopLeft' figzmin];

    [~,cdel] = find(isnan(zbox1));
    zbox1(:,cdel) = [];
    xbox1(:,cdel) = [];
    ybox1(:,cdel) = [];
    fill3(ax2, xbox1,ybox1,zbox1,[.9 .9 .9], 'EdgeColor','k')

    
    
    ZBEDDrawTopRight = ZBEDDrawTop(1+trimNum,colNumberStart: end-trimNum);
    ZBEDDrawTopRight(ZBEDDrawTopRight > zthisH) = zthisH;
    xRight = XSurf(1+trimNum,colNumberStart: end-trimNum);
    yRight = YSurf(1+trimNum,colNumberStart: end-trimNum);
    
    xbox2=[figxmax xthis xRight figxmax];
    ybox2=[figymin figymin yRight figymin];
    zbox2=[figzmin figzmin ZBEDDrawTopRight figzmin];

    [~,cdel] = find(isnan(zbox2));
    zbox2(:,cdel) = [];
    xbox2(:,cdel) = [];
    ybox2(:,cdel) = [];

    fill3(xbox2,ybox2,zbox2,[.7 .7 .7],'EdgeColor','k')
    
    
    xRightBack = XSurf(1+trimNum,1+trimNum: end-trimNum);
    yRightBack = YSurf(end-trimNum,1+trimNum: end-trimNum);
    ZBEDDrawTopRightBack = ZBEDDrawTop(end-trimNum,1+trimNum: end-trimNum);
    
    xbox3=[figxmax figxmin xRightBack figxmax];
    ybox3=[figymax figymax yRightBack figymax];
    zbox3=[figzmin figzmin ZBEDDrawTopRightBack figzmin];

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
      surf(ax1, x3D, xGridNew,zshow,thiscs,'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour
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
      ythis = YSurf(:,ii);
      zshow2 = zthis;
      zshow2(zshow2>zthisH) = NaN;
      
      plot3(ax1, XSurf(:,ii), ythis, zshow2, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);
      
 
      
      % =============== cross section along x ==================
      for j = 1+trimNum : 1+trimNum %size(ZBED_xz,3)
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
          
          xtemp = xGrid_x(1,:);
          [~, coldel] = find(xtemp < xthis);
          
          xGrid_x_new = xGrid_x;
          y3D_new = y3D;
          zshow_new = zshow;
          thiscs_new = thiscs;
          
          xGrid_x_new(:,coldel) = [];
          y3D_new(:,coldel) = [];
          zshow_new(:,coldel) = [];
          thiscs_new(:,coldel) = [];
          
          s = surf(ax1, xGrid_x_new, y3D_new, zshow_new, thiscs_new, 'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour
          % view([0 90]);

          colormap(ax1, mymap);
          caxis manual
          caxis([1 6]);

          hold on;
          ythis = YSurf(j,:);

          zshow2 = zthis;
          zshow2(zshow2>zthisH) = NaN;
          
          XSurfnew = XSurf(j,:);
          [~,coldel] = find(XSurfnew < xthis);
          ythisnew = ythis;
          zshow2new = zshow2;
          XSurfnew(:,coldel) = [];
          ythisnew(:,coldel) = [];
          zshow2new(coldel,:) = [];

          plot3(ax1, XSurfnew, ythisnew, zshow2new, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);
      end
      
      view(-38,20);
      hold off;
      % set(gcf,'GraphicsSmoothing','off');
      figname = sprintf('%s\\xx_%s.png', foldername, num2str(ii-trimNum));
      print('-dpng','-r600',figname);
      
%     figname2 = sprintf('%s\\z_%s.eps', foldername, num2str(ii-trimNum));
%     print(gcf,'-depsc','-painters',figname2); 
%     epsclean(figname2);

      clf;
      
  end


