cd(outputFolder_simulation);   


tic;
colourMatrixR = nan(H,ImageDis,INum);
colourMatrixG = nan(H,ImageDis,INum);
colourMatrixB = nan(H,ImageDis,INum);
sliceTemp = 0;

colourMatrixFacies = nan(H,ImageDis,INum);


for j = 1 : INum

    
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

cd(outputFolder_simulation);
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
trimNum = trimEdge/SurfGridSpace;

figxmin = min(min(XSurf)) + trimEdge;
figxmax = max(max(XSurf)) - trimEdge;
figymin = min(min(YSurf)) + trimEdge;
figymax = max(max(YSurf)) - trimEdge;
figzmin = min(min(min(ZBEDDraw)));
figzmax = max(max(max(ZBEDDraw)));

ZBEDDrawTop = ZBEDDraw(:,:,1);
linesW = 0.1;

zthisH = axisThis_z(1,size(CS_XY, 3));
dis = abs(zrefList - zthisH);
[~, I] = min(dis);
thisHorContLines = ContourLinesList{I};

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
 

 figure;
 fg = figure(1);
      
 set(0,'CurrentFigure',fg);
 set(gcf,'color','w');

    
      
% =============== cross section along x ==================
for j = 1+trimNum : crosssectionBetweenstep : size(ZBED_xz,3)-trimNum
  zthis = ZBED_xz(:,:,j);
  ythis = YSurf(j,1);
  xthis = XSurf(j,1);

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
  
  s = surf(xGrid_x_new, y3D_new, zshow_new, thiscs_new, 'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour

 % axis([figxmin figxmax figymin figymax figzmin figzmax]);
 colormap(mymap);
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

  % plot3(XSurfnew, ythisnew, zshow2new, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);
end


for ii = 1+trimNum : crosssectionBetweenstep : size(ZBED_yz,3)-trimNum

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


    gcf;
 
    surf(x3D, xGridNew,zshow,thiscs,'EdgeColor','none','FaceColor','flat'); % using patch? to fill the predefined colour
    hold on;
    ythis = YSurf(:,ii);
    zshow2 = zthis;
    zshow2(zshow2>zthisH) = NaN;
  
  % plot3(XSurf(:,ii), ythis, zshow2, 'Color', [0.870588 0.721569 0.529412],'Linewidth',linesW);
  

end

 
axis([figxmin figxmax figymin figymax figzmin figzmax]);
daspect([1 1 z_axis_exaggerated]);
grid off;
box on;
% axis off;

view(az,el);

hold off;

xstep = 20;
ystep = 20;
zstep = 5;
fontsizset = 28;
set(gca,'xtick',-40:xstep:40,'GridLineStyle',':'); 
set(gca,'ytick',-40:ystep:40,'GridLineStyle',':');
set(gca,'ztick',-30:30:0,'GridLineStyle',':');

xlabel('x','fontsize',fontsizset,'fontweight','normal','FontAngle','italic','HorizontalAlignment','center');
ylabel('y','fontsize',fontsizset,'fontweight','normal','FontAngle','italic');
zlabel('H','fontsize',fontsizset,'fontweight','normal','FontAngle','italic');
set(gca,'fontsize',fontsizset,'XTickMode','manual','YTickMode','manual','LineWidth',0.5,'TickDir','out');

hZLabel = get(gca,'ZLabel');
set(hZLabel,'rotation',0,'VerticalAlignment','bottom');


m = get(gca,'xtick');
set(gca,'xticklabel',sprintf('%.0f\n',m'),'XMinorTick','off');  
n = get(gca,'ytick');
set(gca,'yticklabel',sprintf('%.0f\n',n'),'YMinorTick','off');

set(gcf,'PaperUnits','centimeters');
outputsize = 30;
set(gcf,'PaperSize',[outputsize outputsize]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperPosition',[0 0 outputsize outputsize]);% [0 0  width height] % unit: cm

 % set(gcf,'GraphicsSmoothing','off');

 figname = sprintf('%s\\fenceDiagram_%s.png', foldername, num2str(1));
 print('-dpng','-r600',figname);
  
    

