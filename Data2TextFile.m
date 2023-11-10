cd(outputFolder_simulation);   


tic;
colourMatrixR = nan(H,ImageDis,INum);
colourMatrixG = nan(H,ImageDis,INum);
colourMatrixB = nan(H,ImageDis,INum);

sliceTemp = 0;
colourMatrixFacies = nan(H,ImageDis,INum);


for j = 1 : INum

    cd(outputFolder_simulation);
   
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
    
    [r1,c1] = find(R==96); % dark gray
    for i = 1 : numel(r1)
        F(r1(i),c1(i)) = 1; 
    end
    
    [r2,c2] = find(R==139); % gray
    for i = 1 : numel(r2)
        F(r2(i),c2(i)) = 2;
    end
    
    [r3,c3] = find(G==190); % orange
    for i = 1 : numel(r3)
        F(r3(i),c3(i)) = 3;
    end
    
    [r4,c4] = find(G==170); % pink
    for i = 1 : numel(r4)
        F(r4(i),c4(i)) = 4;
    end
    
    [r6,c6] = find(G==236); % yellow
    for i = 1 : numel(r6)
        F(r6(i),c6(i)) = 5;
    end

    [r7,c7] = find(G==111); % blue
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



filename = 'Dune_1.txt';
fileID = fopen(filename,'w');
Title = ['Dune ', 'X=', num2str(ImageDis), ' Y=', num2str(INum), ' Z=', num2str(H)]; % PETREL: Properties
properties = 1;
propertyName = 'Facies';


fprintf(fileID,'%s\r\n',Title);
fprintf(fileID,'%d\r\n',properties);
fprintf(fileID,'%s\r\n',propertyName);




for i = H : -1 : 1

        tic;
        for si = 1 : INum
            for k =  1: ImageDis

                fprintf(fileID,'%d\r\n',colourMatrixFacies(i,k,si));

            end
        end
        toc;
        disp(['H ',num2str(i)]);


end
fclose(fileID);

