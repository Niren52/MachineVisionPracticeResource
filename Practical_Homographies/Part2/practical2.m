function r=practical2

%This project explores the geometry of a single camera. The aim is to take several points on
%a plane, and predict where they will appear in the camera image. Based on these observed
%points, we will then try to re-estimate the Euclidean transformation relating the plane and
%the camera. In practical 2b we will use this code to draw a wireframe cube
%on an augmented reality marker.   You should use this
%template for your code and fill in the missing sections marked "TO DO"


%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];
 
%We will assume an object co-ordinate system with the Z-axis pointing upwards and the
%origin in the centre of the plane. There are four known points on the plane, with coordinates
%(in mm):

XCart = [-100 -100  100  100 0 ;...
         -100  100  100 -100 0;...
          0    0    0    0   0 ];

%We will assume that the correct transformation from the plane co-ordinate system to the
%camera co-ordinate system (extrinsic matrix) is:

T = [ 0.9851  -0.0492  0.1619  46.00;...
     -0.1623  -0.5520  0.8181  70.00;...
      0.0490  -0.8324 -0.5518  500.89;...
      0        0       0       1]
  display(size(T));
  
% TO DO  Use the general pin-hole projective camera model discussed in the lectures to estimate 
%where the four points on the plane will appear in the image.  Fill in the
%details of the function "projectiveCamera" - body of function appears below

xImCart = projectiveCamera(K,T,XCart);

% TO DO Add noise to the pixel positions to simulate having to find these points in a noisy
%image. Store the results back in xImCart.  
%The noise should have standard deviation of one pixel in each direction.
xImCart = xImCart + randn(size(xImCart));

%Now we will take the image points and the known positions on the card and try to
%estimate the extrinsic matrix using the algorithm discussed in the lecture. 
%Fill in the details of the function "estimate plane pose" - body of function appears
%below

TEst = estimatePlanePose(xImCart,XCart,K)

%if you have got this correct, it should resemble T above.

%==========================================================================
%==========================================================================

%goal of function is to project points in XCart through projective camera
%defined by intrinsic matrix K and extrinsic matrix T.
function xImCart = projectiveCamera(K,T,XCart)

XHom = [XCart;ones(1,size(XCart,2))];

xExtrinsic = T*XHom;

xCamHom = xExtrinsic(1:3,:);

xImHom = K*xCamHom;

xImCart = xImHom(1:2,:)./repmat(xImHom(3,:),2,1);


%==========================================================================
%==========================================================================

%goal of function is to estimate pose of plane relative to camera
%(extrinsic matrix) given points in image xImCart, points in world XCart
%and intrinsic matrix K.

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


