clear all
close all

A1 = [
    1 0.5 0
    0.5 1 -0.1
    0 -0.1 1
];
center1 = [-1;2;0];

A2 = [
    10 0 0
    0 1 0
    0 0 10
];
center2 = [-1;3;0];

figure(1)
Ellipse_plot(A1, center1)
hold on
Ellipse_plot(A2, center2)

t0 = tic;
[d, x, y] = ...
    ellipsoid_distance( ...
        A1,center1,...
        A2,center2);%,...
%         c1,c2,...
%         epsilon,...
%         max_iter...
%         );
tf = toc(t0);

figure(1)
plot3([x(1),y(1)],[x(2),y(2)],[x(3),y(3)],'r-*')
title(['Distance is ', num2str(d), ' in ', num2str(tf), ' s'])
