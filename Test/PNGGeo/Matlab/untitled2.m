I=imread('./A0.jpg');
I=rgb2gray(I);


I=im2bw(I);
I=~I;
I=imresize(I,0.5);
%  I = bwareaopen(I, 200);
 se = strel('sphere',3);
 I = imdilate(I,se);

  I=imresize(I,5);
%    se = strel('sphere',10);
%     I = imdilate(I,se);
   imshow(~I);
%  I = imerode(I,se);
% % 
I=~I;
  I0=I;
  I=cat(3,I,255.*I0); I=cat(3,I,255.*I0);
  imshow(I);
  imwrite(I,"A.png");