VideoDir = 'exhibition_area';
VideoList = dir(sprintf('%s/*.exr',VideoDir));

nFrames = length(VideoList);
cform = makecform('srgb2lab');
%nFrames = 1;

A = double(exrread(sprintf('%s/%s',VideoDir,VideoList(1).name)));
size_frame = size(A);
A_lab = applycform(A,cform);
ab = double(A_lab(:,:,2:3));
A_reshape = reshape(ab, size_frame(1) * size_frame(2),2);
smooth_lab = ksdensity(A_reshape);
K = size(findpeaks(smooth_lab),1);
%K = 30;

image_sequence = zeros([size_frame,nFrames]);

context.K = K;
context.initialize = true;
context.initial_iters = 10;
context.iters = 2;
context.sigma = 1/200 * max(size_frame(1),size_frame(2));

videoWriter = VideoWriter('/Users/ntaheria/Desktop/output33.avi');
open(videoWriter);

for i = 100 : nFrames
    fprintf('Frame %d\n', i)
    frame = exrread(sprintf('%s/%s',VideoDir,VideoList(i).name));
    frame = double(frame);
    %max_frame = max(max(max(frame)));
    %frame = frame/max_frame;
    [frame_tonemapped, context] = tone_mapping(frame,context);
    image_sequence(:,:,:,i) = frame_tonemapped;
    video_frame = min(max(frame_tonemapped, 0), 1);
    writeVideo(videoWriter, video_frame);
end
%imshow(frame_tonemapped);
close(videoWriter);
imshow([frame, frame_tonemapped])