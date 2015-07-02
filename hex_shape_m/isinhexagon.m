function [ inhexagon ] = isinhexagon(x,y,d)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tol=1.001;

%tol=1.05;
 
war1 = (abs(y) <= d*tol*sqrt(3)/2);
war2 = (abs(sqrt(3)/2*x+1/2*y)) <= tol*d*sqrt(3)/2;
war3 = (abs(sqrt(3)/2*x-1/2*y)) <= tol*d*sqrt(3)/2;

inhexagon=war1&war2&war3;

