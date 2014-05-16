clear;clc;
% make sure we add the correct folders even if this file is
% not called from the current folder
fileName = mfilename();
filePath = mfilename('fullpath');
filePath = filePath(1:end-size(fileName, 2));

% Add folders to current path
path(genpath([filePath 'Files']), path);

fprintf('Automatic Text Detection: Starting DEMO...\n');

% Open main GUI
TextDetectionDemo;

fprintf('Automatic Text Detection: Ready.\n\n');

% clear variables
clearvars fileName filePath