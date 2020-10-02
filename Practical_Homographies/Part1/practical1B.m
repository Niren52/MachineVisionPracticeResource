function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');
figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2
H1to2 = calcBestHomeography(pts1, pts2);
[w1,h1,~] = size(im1);
[w2,h2,~] = size(im2);
[w3,h3,~] = size(im3);
[X,~] = meshgrid(1:h1,1:w1);
idx_i = repmat(1:w1,[1,h1]);
idx_j = reshape(X,[1,w1*h1]);
tr_pos = H1to2*[idx_j;idx_i;ones(1,size(im1,1)*size(im1,2))];
tr_pos(1:2,:) = floor(tr_pos(1:2,:)./repmat(tr_pos(3,:),[2,1]));
%Routine:
for j=1:h1
    for i=1:w1
        new_x = tr_pos(2,i+(j-1)*w1);
        new_y = tr_pos(1,i+(j-1)*w1);
        if(new_x > 0 && new_x < w2 ...
        && new_y > 0 && new_y < h2)
            im1(i,j,:) = im2(new_x,new_y,:);
        end
    end
end

%****TO DO****
%repeat the above process mapping image 3 to image 1.
H1to3 = calcBestHomeography(pts1b, pts3);
tr_pos = H1to3*[idx_j;idx_i;ones(1,size(im1,1)*size(im1,2))];
tr_pos(1:2,:) = floor(tr_pos(1:2,:)./repmat(tr_pos(3,:),[2,1]));

%Routine:
for j=1:h1
    for i=1:w1
        new_x = tr_pos(2,i+(j-1)*w1);
        new_y = tr_pos(1,i+(j-1)*w1);
        if(new_x > 0 && new_x < w3 ...
        && new_y > 0 && new_y < h3)
            im1(i,j,:) = im3(new_x,new_y,:);
        end
    end
end
%Display figure
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;


%TAKEN FROM PRACTICAL 1
function H = calcBestHomeography(pts1Cart, pts2Cart)

pts1Cart = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Cart = [pts2Cart; ones(1,size(pts2Cart,2))];

matrixA = zeros(2*size(pts1Cart,2),9);
for i = 1:size(pts1Cart,2)
    a = pts1Cart(1,i);
    b = pts1Cart(2,i);
    c = pts2Cart(1,i);
    d = pts2Cart(2,i);
    matrixA(2*i-1,:) = [0,0,0,-a,-b,-1,d*a,d*b,d];
    matrixA(2*i,:) = [a,b,1,0,0,0,-c*a,-c*b,-c];
end

h = solveAXEqualsZero(matrixA); 

H = reshape(h,[3,3])';

%TAKEN FROM PRACTICAL 1
function x = solveAXEqualsZero(matrixA)
[~,~,V] = svd(matrixA);
x = V(:,end);

