% Code for creating the pointcloud format file from the 3D point-clouds 
% genertated from image3D_integrate.m

clc;
clear;

close all;

fileSave = 'parkingLot_full_L1.ply';
load('resultsCombine');
height_threshol_above = 7; %  +7m above ground
height_threshold_below = -1; %  -1m above ground


idx=find(pointsTotalFinal(:,3)>height_threshol_above);
pointsTotalFinal(:,idx) = [];
ampsTotalFinal(:,idx) = [];

idx=find(pointsTotalFinal(:,3)<height_threshold_below);
pointsTotalFinal(:,idx) = [];
ampsTotalFinal(:,idx) = [];

pointCloud= pointCloud(pointsTotalFinal);
normals = pcnormals(pointCloud,12);
pointCloud.Normal=normals;


pcwrite(pointCloud,fileSave);

