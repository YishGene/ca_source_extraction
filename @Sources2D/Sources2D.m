classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data. 
    % Author: Pengcheng Zhou, zhoupc1988@gmail.com with modifications from
    % Eftychios Pnevmatikakis
    
    %% properties 
    properties
        A;          % spatial components of neurons 
        C;          % temporal components of neurons 
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        S;          % spike counts 
        Coor;       % neuron contours         
        Df;         % background for each component to normalize the filtered raw data  
        C_df;       % temporal components of neurons and background normalized by Df        
        S_df;       % spike counts of neurons normalized by Df
        options;    % options for model fitting 
        P;          % some estimated parameters 
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = Sources2D(varargin)
            obj.options = CNMFSetParms(); 
            obj.P = struct('p', 2); 
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:}); 
            end
        end 
        
        %% update parameters
        function updateParams(obj, varargin)
                obj.options = CNMFSetParms(obj.options, varargin{:});          
        end
        
        %% data preprocessing
        function Y = preprocess(obj,Y,p)
            [obj.P,Y] = preprocess_data(Y,p,obj.options);
        end
        
        %% fast initialization
        function [center, res] = initComponents(obj, Y, K, tau)
            [obj.A, obj.C, obj.b, obj.f, center, res] = initialize_components(Y, K, tau, obj.options);
        end
        
        %% add new components post-initialization
        % prev_res is a residual left from previous calculations (i.e.
        % minus the components already accounted for. Saves the need to
        % recalculate the residual
        function [center, res] = addComponents(obj, prev_res, K, tau)
            [newA, newC, ~, ~, center, res] = initialize_components(prev_res, K, tau, obj.options);
            obj.A = [obj.A, newA];
            obj.C = [obj.C; newC];
        end
        
        
        
        %% manual refinement
        function center = refineComponents(Y,obj,center,img,sx)
            oldDoSeedROI = obj.options.doSeedROI;
            obj.options.doSeedROI = false;
                
            [obj.A,obj.C,center] = manually_refine_components(Y,obj.A,obj.C,center,img,sx,obj.options);
            
            obj.options.doSeedROI = oldDoSeedROI;
        end
        
        %% update spatial components
        function updateSpatial(obj, Y)
            [obj.A, obj.b, obj.C] = update_spatial_components(Y, ...
                obj.C, obj.f, obj.A, obj.P, obj.options);
        end
        
        %% update temporal components
        function updateTemporal(obj, Y, varargin)
            % updateTemporal(obj, Y, ..., parameter, value);
            %
            % parameter values can be of the following: 
            % 'A', 'b', 'C', 'f', 'P', 'options'
            % 
            % default values are taken from obj
            %
            % example: obj.updateTemporal(Y, 'A', Aprevious, 'C', Cprevious)
            
            
            parser = inputParser;
            parser.addParameter('A', obj.A);
            parser.addParameter('b', obj.b);
            parser.addParameter('C', obj.C);
            parser.addParameter('f', obj.f);
            parser.addParameter('P', obj.P);
            parser.addParameter('options', obj.options);
            
            parser.parse(varargin{:});
            
            A_ = parser.Results.A;
            b_ = parser.Results.b;
            C_ = parser.Results.C;
            f_ = parser.Results.f;
            P_ = parser.Results.P;
            options_ = parser.Results.options;
            
            [obj.C, obj.f, obj.P, obj.S] = update_temporal_components(...
                Y, A_, b_, C_, f_, P_, options_);
        end
                       
        %% merge found components
        function [nr, merged_ROIs] = merge(obj, Y, varargin)
            % inputs:
            %   obj:    
            %   Y  : the imaging matrix
            %   force_merge_groups (optional): a cell array of groups to
            %   use when force-merging ROIs
            [obj.A, obj.C, nr, merged_ROIs, obj.P, obj.S] = merge_components(...
                Y,obj.A, obj.b, obj.C, obj.f, obj.P,obj.S, obj.options, varargin{:});
        end
        
        %% 
        function  selectandmergeROIs(obj, Y, meanstack)
            if ~exist('meanstack', 'var'), 
                meanstack = [];
            end
            center = com(obj.A,obj.options.d1,obj.options.d2);
            hsc = scatter(center(:,2), center(:,1));
            
            rect = getrect; % get a rectangle that contains the centers of the suspect ROIs
            tf = arrayfun(@(x,y)(x>rect(1) && x< rect(1)+rect(3) && y>rect(2) && y < rect(2)+rect(4)), center(:,2), center(:,1));
            findtf = find(tf);
            if isempty(tf), 
                return;
            end
            obj.plotComponents(Y, meanstack, findtf); 
            uiwait(gcf)
            
            answer = questdlg('Merge components?', 'Merge question');
            switch answer, 
                case 'Yes'
                    obj.merge(Y, {findtf});
            end
            delete(hsc)
        end
        %% compute the residual
        function [Y_res] = residual(obj, Yr)
            Y_res = Yr - obj.A*obj.C - obj.b*obj.f;
        end
        
        %% take the snapshot of current results
        function [A, C, b, f, P] = snapshot(obj)
            A = obj.A;
            C = obj.C;
            b = obj.b;
            f = obj.f;
            P = obj.P;
        end
        
        %% extract DF/F signal after performing NMF
        function [C_df, Df, S_df] = extractDF_F(obj, Y, i)
            if ~exist('i', 'var')
                i = size(obj.A, 2) + 1;
            end
            
            [obj.C_df, obj.Df, obj.S_df] = extract_DF_F(Y, [obj.A, obj.b],...
                [obj.C; obj.f], obj.S, i);
            
            C_df =  obj.C_df;
            Df = obj.Df;
            S_df = obj.S_df;
            
        end
        
        %% order_ROIs 
        function [srt] = orderROIs(obj)
            [obj.A, obj.C, obj.S, obj.P, srt] = order_ROIs(obj.A, obj.C,...
                obj.S, obj.P); 
        end 
        
        %% view contours 
        function [json_file] = viewContours(obj, Cn, contour_threshold, display)
            if ~exist('Cn', 'var') || isempty(Cn)
                Cn = reshape(obj.P.sn, obj.options.d1, obj.options.d2); 
            end
            if ~exist('contour_threshold', 'var') || isempty(contour_threshold) 
                contour_threshold = obj.options.cont_threshold; 
            end
            if ~exist('display', 'var')  || isempty(display)
                display = 0; 
            end
            [obj.Coor, json_file] = plot_contours(obj.A, Cn, ...
                contour_threshold, display); 
        end 
        
        %% plot components 
        function plotComponents(obj, Y, Cn, indROIs)
            if ~exist('Cn', 'var')
                Cn = []; 
            end
            if ~exist('indROIs', 'var'), 
                indROIs = [];
            end
            if ~isempty(obj.Df), % save time by avoiding recalculation of Df
                view_components(Y, obj.A, obj.C, obj.b, obj.f, Cn, obj.options, indROIs, obj.Df); 
            else
                view_components(Y, obj.A, obj.C, obj.b, obj.f, Cn, obj.options, indROIs); 
            end
        end 
        
        %% plot components GUI
        function plotComponentsGUI(obj, Y, Cn)
            if ~exist('Cn', 'var')
                Cn = []; 
            end
            plot_components_GUI(Y,obj.A,obj.C,obj.b,obj.f,Cn,obj.options)
        end 
        
        %% make movie 
        function makePatchVideo(obj, Y) 
            make_patch_video(obj.A, obj.C, obj.b, obj.f, Y, obj.Coor,...
                obj.options); 
        end 
    end
    
end