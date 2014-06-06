function [ i,j ] = getIJ( index,w,height)
%GETIJ Return the coordinates of the patch given its index in the matrix
%obtained by im2col
%   Detailed explanation goes here
i=w+mod(index-1,height-2*w)+1;
j=1+w+floor((index-1)/(height-2*w));
end
