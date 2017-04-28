VideoDir = 'bridge_2';
VideoList = dir(sprintf('%s/*.exr',VideoDir));

%nFrames = length(VideoList);
nFrames = 20;

A = double(exrread(sprintf('%s/%s',VideoDir,VideoList(1).name)));
size_frame = size(A);
A_reshape = reshape(A, size_frame(1) * size_frame(2),3);
max_A = max(max(max(A_reshape)));
A_reshape = A_reshape/max_A;
hsv = rgb2hsv(A_reshape);
smooth_hue_hist = ksdensity(hsv(:,1));
K = size(findpeaks(smooth_hue_hist),2);
image_sequence = zeros([size_frame,nFrames]);

context.K = K;
context.initialize = true;
context.initial_iters = 6;
context.iters = 1;
context.sigma = 1/200 * max(size_frame(1),size_frame(2));

videoWriter = VideoWriter('/Users/ntaheria/Desktop/output14.avi');
open(videoWriter);

for i = 1 : nFrames
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

close(videoWriter);