function inflectionPointsXYCol = InflectionPoints(lineX,lineY, pns)

InflectionsNum = 1;
for i = 3: pns-3 %2 : pns - 2
  % First take 2 points on the current contour :
  r0 = [lineX(i); lineY(i)];
  r1 = [lineX(i+1); lineY(i+1)];
  % let them define a stratight line on the form: 
  % dot(r,e_n) - 1 = 0
  e_n = [cos(pi/2) sin(pi/2);-sin(pi/2) cos(pi/2)]*(r1-r0)/norm(r1-r0);
  l = dot(r0,e_n);
  % the points on the contour on either side of the line-segment:
  r_p = [lineX(i+2); lineY(i+2)];
  r_m = [lineX(i-1); lineY(i-1)];
  
  r_p1 = [lineX(i+3); lineY(i+3)];
  r_m1 = [lineX(i-2); lineY(i-2)];

  % lengths along e_n to points r_p and r_m
  l_p = dot(r_p,e_n);
  l_m = dot(r_m,e_n);
  
  l_p1 = dot(r_p1,e_n);
  l_m1 = dot(r_m1,e_n);
  
  % if they are on either side of the line-segment then we have
  % an inflection point
  if (l_p - l)*(l_m - l) <0 && (l_p1 - l)*(l_m1 - l) < 0
    % take the mid-point of that line-segment as the
    % inflection-point, change this according to taste...
    inflectionPointsXY(InflectionsNum,:) = (r0+r1)'/2;
    inflectionPointsCol(InflectionsNum,1) = i;
    InflectionsNum = InflectionsNum + 1;
  end
end
if exist('inflectionPointsXY','var')
    inflectionPointsXYCol = [inflectionPointsXY,inflectionPointsCol];
else
    inflectionPointsXYCol = [];
end
