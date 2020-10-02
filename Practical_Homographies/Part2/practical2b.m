function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"
close all
clear all
%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.
%Routine function below
TEst = estimatePlanePose(xImCart,XCart,K);

%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];
lineProjectionA = [K, zeros(3,1)]*TEst*[XWireFrameCart;ones(1,size(XWireFrameCart,2))];
%CALCULATE PROJECTED POSITIONS
lineProjectionB = lineProjectionA(1:3,:)./repmat(lineProjectionA(3,:),3,1);

display(lineProjectionB);
display(lineProjectionB(1,1));

%JOIN PROJECTED POSITIONS TOGETHER FOR EACH POSITION
plot([lineProjectionB(1,1),lineProjectionB(1,2)],[lineProjectionB(2,1),lineProjectionB(2,2)],'g-')
plot([lineProjectionB(1,1),lineProjectionB(1,4)],[lineProjectionB(2,1),lineProjectionB(2,4)],'g-')
plot([lineProjectionB(1,1),lineProjectionB(1,5)],[lineProjectionB(2,1),lineProjectionB(2,5)],'g-')
plot([lineProjectionB(1,2),lineProjectionB(1,6)],[lineProjectionB(2,2),lineProjectionB(2,6)],'g-')
plot([lineProjectionB(1,2),lineProjectionB(1,3)],[lineProjectionB(2,2),lineProjectionB(2,3)],'g-')
plot([lineProjectionB(1,3),lineProjectionB(1,4)],[lineProjectionB(2,3),lineProjectionB(2,4)],'g-')
plot([lineProjectionB(1,3),lineProjectionB(1,7)],[lineProjectionB(2,3),lineProjectionB(2,7)],'g-')
plot([lineProjectionB(1,4),lineProjectionB(1,8)],[lineProjectionB(2,4),lineProjectionB(2,8)],'g-')
plot([lineProjectionB(1,5),lineProjectionB(1,8)],[lineProjectionB(2,5),lineProjectionB(2,8)],'g-')
plot([lineProjectionB(1,5),lineProjectionB(1,6)],[lineProjectionB(2,5),lineProjectionB(2,6)],'g-')
plot([lineProjectionB(1,6),lineProjectionB(1,7)],[lineProjectionB(2,6),lineProjectionB(2,7)],'g-')
plot([lineProjectionB(1,7),lineProjectionB(1,8)],[lineProjectionB(2,7),lineProjectionB(2,8)],'g-')

                
%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?


%==========================================================================
%==========================================================================

%goal of function is to estimate pose of plane relative to camera
%(extrinsic matrix) given points in image xImCart, points in world XCart
%and intrinsic matrix K.

%FROM PART 2A

function T = estimatePlanePose(xImCart,XCart,K)
xImHom = [xImCart;ones(1,size(xImCart,2))];

xCamHom = K^(-1)*xImHom;

H = calcBestHomography(XCart,xCamHom);

Phi = H(:,1:2);
[U,~,V] = svd(Phi);
L = [1 0; 0 1; 0 0];
R = U*L*V';
R = [R, cross(R(:,1),R(:,2))];

if(det(R)<0)
    R(:,3) = -R(:,3);
end

l = sum(sum(H(:,1:2)./R(:,1:2)))/6;
t = H(:,3)/l;

if(t(3)<0)
    t=-t;
    R(:,1:2) = -R(:,1:2);
end

T  = [R t;0 0 0 1];


%TAKEN FROM PART 1
function H = calcBestHomography(pts1Cart, pts2Cart)

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

function x = solveAXEqualsZero(matrixA)
[~,~,V] = svd(matrixA);
x = V(:,end);


