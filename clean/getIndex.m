function [index]=getIndex(i,j,w,height)
%GETINDEX Return the index of the patch given its coordinates i and j.
% i and j should satisfy i>=1+w and i<=m-w, j>=1+w;
index=(j-w-1)*(height-2*w)+i-w;
end

