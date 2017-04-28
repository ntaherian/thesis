VideoDir = 'bridge_2';
VideoList = dir(sprintf('%s/*.exr',VideoDir));

nFrames = length(VideoList);
cform = makecform('srgb2lab');
nFrames = 10;

% A = double(exrread(sprintf('%s/%s',VideoDir,VideoList(1).name)));
% size_frame = size(A);
% A_lab = applycform(A,cform);
% ab = double(A_lab(:,:,2:3));
% A_reshape = reshape(ab, size_frame(1) * size_frame(2),2);
% smooth_lab = ksdensity(A_reshape);
% context.K = size(findpeaks(smooth_lab),1);

% A_reshape = reshape(A, size_frame(1) * size_frame(2),3);
% max_A = max(A_reshape(:));
% A_reshape = A_reshape/max_A;
% hsv = rgb2hsv(A_reshape);
% smooth_hue_hist = ksdensity(hsv(:,1));
% context.K = size(findpeaks(smooth_hue_hist),2);

image_sequence = zeros([size_frame,nFrames]);

context.initialize = true;
context.initial_iters = 10;
context.iters = 1;
context.sigma = 1/200 * max(size_frame(1),size_frame(2));

videoWriter = VideoWriter('/Users/ntaheria/Desktop/output61.avi');
open(videoWriter);
distances = [];

for i = 1 : nFrames
    fprintf('Frame %d\n', i)
    frame = exrread(sprintf('%s/%s',VideoDir,VideoList(i).name));
    frame = double(frame);
    if (i > 1)
        distance = scene_change(frame,prev_frame)
        distances = [distances distance];
        if (distance > 1.6)
            context.initialize = true; 
        end
    end
    prev_frame = frame;

    %max_frame = max(frame(:));
    %frame = frame/max_frame;
    
    [frame_tonemapped, context] = tone_mapping(frame,context);
    image_sequence(:,:,:,i) = frame_tonemapped;
    video_frame = min(max(frame_tonemapped, 0), 1);
    writeVideo(videoWriter, video_frame);
end
%imshow(frame_tonemapped);
close(videoWriter);
%imshow([frame, frame_tonemapped])
plot(distances)