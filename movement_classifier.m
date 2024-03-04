% Movement classification:
clear all; close all; clc;

% Reading folder with videos to train the classifier
datasetFolder = 'C:\Users\Atti\Desktop\Crina\all';
classNames = {'walk','jump', 'run', 'altc'};
videoFiles = dir(fullfile(datasetFolder, '*.avi'));

features = [];
labels = [] ;

% defining the video player
videoPlayer = vision.VideoPlayer();

% iterating trough each video
for i = 1:length(videoFiles)  
    % defining the detectors for foreground extraction
    [detector, detector1, blob, blob1, blobFeet]  = detector_init();
    
    % reading the video
    filename = fullfile(datasetFolder, videoFiles(i).name);
    videoSource = VideoReader(filename);

    aux_f = [];

    % iterating trough each frame
    while hasFrame(videoSource)
        clear s;
        % reading the frame and generating mask
        frame  = readFrame(videoSource);
        fgMask = detector(frame);
        out    = frame;
        
        % gettting the features
        if ~isempty(blob(fgMask))
            
            centroids = uint8(blob(fgMask)); % the hips centroid
            
    
            feet = uint8(blob1(fgMask)); % box around the feet
            feet =[feet(1)-5, feet(2)+feet(4)/2, feet(3)+5, feet(4)/2]; 
           
            s = regionprops(fgMask,'Extrema'); % geting the extremas from the blops
            k = s.Extrema;
    
            % selectring the extremas from inside the box around the feet
            for k = 1:length(s)
                blob_aux(:,:,k) = s(k).Extrema;
                
                cond1 = ((blob_aux(4,1,k) >= feet(1)) && (blob_aux(4,1,k) <= feet(1) + feet(3)));
                cond2 = ((blob_aux(4,2,k) >= feet(2)) && (blob_aux(4,1,k) <= feet(2) + feet(4)));
                if cond1 == 1
                x = [ blob_aux(4,1,k) blob_aux(6,1,k)];
                y = [ blob_aux(4,2,k) blob_aux(6,2,k)];
                end
            end
        
            
            % constructing the auxiliar features array which contains the x,y
            % coordinates of the hip and two extremas for each frame
            aux_f = [aux_f; centroids(1,:), x, y];
       
            % constructing the labels
            label = 'altc';
            [~, name, ~] = fileparts(videoFiles(i).name);
            for j = 1:length(classNames)
                cond = contains(lower(name), lower(classNames{j}));
                    if cond == 1
                        % display(contains(lower(name), lower(classNames{j})))
                        label = classNames{j};
                        %display("shit")
                    end
            end
            % padding if needed
            if strcmp("run", label)
                label = '0run';
            end
            labels = [labels; label];
            
            % iinserting the markings on the frame  (visual purpose only)
            out = insertMarker(out,centroids(1,:),"x-mark","Size",5);
            out = insertShape(out,"rectangle",feet);
            
            out = insertMarker(out,[x(1),y(1)],"x-mark","Size",5);
            out = insertMarker(out,[x(2),y(2)],"x-mark","Size",5);
            out = insertShape(out,'line',[x(1),y(1) centroids(1,:)],'LineWidth',3);
            out = insertShape(out,'line',[x(2),y(2) centroids(1,:)],'LineWidth',3);
        
        end
        
        % applying the mask on the video (visual purpose only)
        out(~fgMask) = 0;
        % visualizing the features on the vide
        videoPlayer(fgMask);
        pause(0.01);
    
    end
    %  constructing the features matrix
    features = [features; aux_f];

end

% constructing the final features matrix from the distance between the
% positions between each consecutive frame
x_dist_h =calc_dist(features(:,1)) 
y_dist_h =calc_dist(features(:,2)) 
x_dist_r =calc_dist(features(:,3)) 
y_dist_r =calc_dist(features(:,5))
x_dist_l =calc_dist(features(:,4)) 
y_dist_l =calc_dist(features(:,6))

FEATURES = [x_dist_h, y_dist_h, x_dist_r, y_dist_r,x_dist_l, y_dist_l];

%Storing the obtained values so that we can redefine the classifier, if needed, without obtaining the features again 
FEATURES_CONST = FEATURES;
labels_const = labels;

%% Defining the classidfier
classifier = fitcecoc(FEATURES_CONST, labels_const,'CategoricalPredictors','all');

%% Testing the classifier
clearvars -except classifier FEATURES_CONST labels_const

% The new video to be tested:
newVideoFile = 'C:\Users\Atti\Desktop\Crina\all\denis_jump.avi';
filename = newVideoFile;
videoSource = VideoReader(filename);
video = videoSource;

% defining the detectors for foreground extraction
[detector, detector1, blob, blob1, blobFeet]  = detector_init();
videoPlayer = vision.VideoPlayer();

aux_f = [];
features = [];
% iterating trough each frame
while hasFrame(videoSource)
    
    % reading the frame and generating mask
    frame  = readFrame(videoSource);
    fgMask = detector(frame);
    out    = frame;
    
    % gettting the features
    if ~isempty(blob(fgMask))
        
        centroids = uint8(blob(fgMask)); % the hips centroid
        

        feet = uint8(blob1(fgMask)); % box around the feet
        feet =[feet(1)-5, feet(2)+feet(4)/2, feet(3)+5, feet(4)/2]; 
       
        s = regionprops(fgMask,'Extrema'); % geting the extremas from the blops
        k = s.Extrema;

        % selectring the extremas from inside the box around the feet
        for k = 1:length(s)
            blob_aux(:,:,k) = s(k).Extrema;
            
            cond1 = ((blob_aux(4,1,k) >= feet(1)) && (blob_aux(4,1,k) <= feet(1) + feet(3)));
            cond2 = ((blob_aux(4,2,k) >= feet(2)) && (blob_aux(4,1,k) <= feet(2) + feet(4)));
            if cond1 == 1
            x = [ blob_aux(4,1,k) blob_aux(6,1,k)];
            y = [ blob_aux(4,2,k) blob_aux(6,2,k)];
            end
        end
    
        
        % constructing the auxiliar features array which contains the x,y
        % coordinates of the hip and two extremas for each frame
        aux_f = [aux_f; centroids(1,:), x, y];
   
            
        
        % iinserting the markings on the frame  (visual purpose only)
        out = insertMarker(out,centroids(1,:),"x-mark","Size",5);
        out = insertShape(out,"rectangle",feet);
        
        out = insertMarker(out,[x(1),y(1)],"x-mark","Size",5);
        out = insertMarker(out,[x(2),y(2)],"x-mark","Size",5);
        out = insertShape(out,'line',[x(1),y(1) centroids(1,:)],'LineWidth',3);
        out = insertShape(out,'line',[x(2),y(2) centroids(1,:)],'LineWidth',3);
    
    end
    
    % applying the mask on the video (visual purpose only)
    out(~fgMask) = 0;
    % visualizing the features on the vide
    videoPlayer(out);
    pause(0.01);

end
%  constructing the features matrix
features = [features; aux_f];

% constructing the final features matrix from the distance between the
% positions between each consecutive frame
x_dist_h =calc_dist(features(:,1)) 
y_dist_h =calc_dist(features(:,2)) 
x_dist_r =calc_dist(features(:,3)) 
y_dist_r =calc_dist(features(:,5))
x_dist_l =calc_dist(features(:,4)) 
y_dist_l =calc_dist(features(:,6))

FEATURES = [x_dist_h, y_dist_h, x_dist_r, y_dist_r,x_dist_l, y_dist_l];

% Introducing the features in the classifier so we obtain a model
predictedClass = predict(classifier, FEATURES);
display(predictedClass(1,:))

%% Functions declared:
% the following function calculates the eucladioan distance between the
% consecutive positions
function out = calc_dist(in)
    obspoints =[in(:) circshift(in(:),1)];
    
    d = pdist(obspoints(2:end,:));
    d_s = squareform(d);

    out = [];
    for m = 1:length(in)-2
        out(m) = d_s(m,m+1);
    end
    out =[out 0 0]'; 

end

% the following function initializes the foreground detectors and
% blobanalyzers
function [detector, detector1, blob, blob1, blobFeet]  = detector_init()
    detector = vision.ForegroundDetector('NumTrainingFrames', 10,'InitialVariance', 30*30);
    detector1 = vision.ForegroundDetector('NumTrainingFrames', 10,'InitialVariance', 30*30);
    blob = vision.BlobAnalysis('CentroidOutputPort', true, 'AreaOutputPort', false, 'BoundingBoxOutputPort', false,'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 400);
    blob1 = vision.BlobAnalysis('CentroidOutputPort', false, 'AreaOutputPort', false, 'BoundingBoxOutputPort', true,'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 400);
    blobFeet = vision.BlobAnalysis('CentroidOutputPort', true, 'AreaOutputPort', false,'MajorAxisLengthOutputPort', true, 'BoundingBoxOutputPort', false,'MinimumBlobAreaSource', 'Property', 'MinimumBlobArea', 200);
end
