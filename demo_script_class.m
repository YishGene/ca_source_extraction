clear;
% same demo as demo_script.m but using the class @Sources2D
%% load file

addpath(genpath('utilities'));
             
nam = 'demoMovie.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double
Y = Y - min(Y(:));                                  % make data non-negative

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% load saved registered data
load('FF_000-002.mat_stackdata', '-mat')
%%
Y = double(allstackreg{1}{1});
Y = Y - min(Y(:));                                  % make data non-negative

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;        


mstackreg = mean(allstackreg{1}{1},3);
mstackregred = mean(allstackregorig{1}{1},3);

roiSeeds = [];

%% find seeds from red channel
clear UpdateCentersOnFigure % clear persistent variables in the function

qr = quantile(mstackregred(1:10:end), [.05 .99]);
qg = quantile(mstackreg(1:10:end), [.05 .99]);
meanstack = cat(3, (mstackregred-qr(1))/diff(qr), (mstackreg-qg(1))/diff(qg), zeros(size(mstackreg), 'like', mstackreg));
figure
him = imshow(meanstack); 
hax = get(him, 'parent');


fred = medfilt2(mstackregred, [5 5]);
figure
himfred = imshow(fred, []);
haxfred = get(himfred, 'parent');
UpdateCentersOnFigure( 1000, fred, hax ); % plot the centers on the image and prep the image for updating. 

addlistener(gca, 'CLim', 'PostSet', @(src,event)UpdateCentersOnFigure(event.AffectedObject.CLim(1)))

imcontrast(haxfred)
msgbox('Set image threshold using the histogram', 'Image threshold', 'modal');


%% Get Centers  
roiSeeds = UpdateCentersOnFigure([]);

%% Set parameters

K = size(roiSeeds,1);                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

obj = Sources2D;
updateParams(obj,...            
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'seedROI', roiSeeds(:, [2 1]),...
    'doSeedROI', true ...
    );

%% Data pre-processing

Y = preprocess(obj,Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[center, res] = initComponents(obj, Y, K, tau);

% display centers of found components
Cn =  reshape(obj.P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% manually refine components (optional)
img = meanstack;
updateParams(obj, 'doSeedROI', false)
refine_components = true;  % flag for manual refinement
if refine_components
    [center] = refineComponents(Y,obj,center,img,tau);
end

%% add more spatial components from residual


    
%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
updateSpatial(obj, Yr);

%% update temporal components
updateTemporal(obj, Yr);

%% merge found components
Apr = obj.A;    % store non-merged components
Cpr = obj.C;
[K_m, merged_ROIs] = merge(obj, Yr);
display_merging = 1; % flag for displaying merging example
if display_merging
    i = 1; randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(Apr(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(obj.A(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(Cpr(merged_ROIs{i},:),[],2))\Cpr(merged_ROIs{i},:))'); 
            hold all; plot(1:T,obj.C(K_m-length(merged_ROIs)+i,:)/max(obj.C(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
updateSpatial(obj, Yr);
updateTemporal(obj, Yr);

%% do some plotting
[srt] = orderROIs(obj);     % order components
K_m = size(obj.C,1);
[C_df, ~, S_df] = extractDF_F(obj, Yr, K_m+1); % extract DF/F values.

contour_threshold = 0.85;   % amount of energy used for each component to construct contour plot
figure;
[json_file] = viewContours(obj, Cn, contour_threshold, 1);
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

plotComponentsGUI(obj, Yr, Cn);     % display all components
pause;

%% plot the raw dF/F for all components separated. 
Y_r = (obj.A'*Yr- (obj.A'*obj.A)*obj.C - (obj.A'*full(obj.b))*obj.f) + obj.C;
plotSpacing(1:14273, bsxfun(@rdivide,Y_r', Df(1:end-1)'), 2, 'k')


%% make movie
makePatchVideo(obj, Yr) 