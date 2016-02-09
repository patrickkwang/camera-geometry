classdef CameraMatrix % < hgsetget
	% class for handling camera transformations
	
	properties
		% default to the identity transformation
		K = eye(3);
		R = eye(3);
		c = zeros(3,1);
		% t = R*c
		imageSize
	end
	
	properties (Dependent)
		x0
		y0
		f
		
		t % = -Rc
		P
	end
	
	methods
		% Allow for string, value pairs
		function obj = CameraMatrix(varargin)
			obj = prtUtilAssignStringValuePairs(obj,varargin{:});
		end
		
% 		function varargout = imageToWorld(obj, varargin)
% 			if nargin==2
% 				Vi = [varargin{1}, ones(size(varargin{1},1),1)]';
% 			elseif nargin==4
% 				Vi = [varargin{1}, varargin{2}, 1]';
% 			end
% 			
% 			% assuming ground plane is z=0
% 			Vw = obj.P(:,[1,2,4])\Vi;
% 			Xw = Vw(1,:)./Vw(3,:);
% 			Yw = Vw(2,:)./Vw(3,:);
% 			
% 			if nargout==1
% 				varargout{1} = [Xw,Yw,0];
% 			elseif nargout==3
% 				varargout{1} = Xw;
% 				varargout{2} = Yw;
% 				varargout{3} = 0;
% 			end
% 		end
		
		function varargout = imageToWorld(obj, varargin)
			if nargin==2
				Vi = [varargin{1}, ones(size(varargin{1},1),1)]';
				Zw = zeros(1,size(varargin{1},1));
			elseif nargin==3
				Vi = [varargin{1}, ones(size(varargin{1},1),1)]';
        if isscalar(varargin{2})
          Zw = varargin{2}*ones(size(varargin{1},1),1)';
        else
          Zw = varargin{2};
        end
			elseif nargin==4
				Vi = [varargin{1}, varargin{2}, 1]';
				Zw = 0;
			elseif nargin==5
				Vi = [varargin{1}, varargin{2}, 1]';
				Zw = varargin{3};
			end
			
			q = obj.R.R'/obj.K*Vi;
			Xc = bsxfun(@times,bsxfun(@rdivide,q,q(3,:)),Zw-obj.c(3));
			Vw = bsxfun(@plus,Xc,obj.c);
			Xw = Vw(1,:);
			Yw = Vw(2,:);
			
			if nargout<2
				varargout{1} = [Xw;Yw;Zw]';
			elseif nargout==3
				varargout{1} = Xw';
				varargout{2} = Yw';
				varargout{3} = Zw';
			end
		end
		
		function varargout = worldToImage(obj, varargin)
			% input should be n-by-d
			% output is n-by-(d-1)
			
			if nargin==2
				Vw = [varargin{1}, ones(size(varargin{1},1),1)]';
			elseif nargin==4
				Vw = [varargin{1}, varargin{2}, varargin{3}, 1]';
			end
			
			Vi = obj.P*Vw;
			Xi = Vi(1,:)./Vi(3,:);
			Yi = Vi(2,:)./Vi(3,:);
			
			if nargout<2
				varargout{1} = [Xi;Yi]';
			elseif nargout==2
				varargout{1} = Xi';
				varargout{2} = Yi';
			end
		end
		
		function plot(obj,varargin)
			colors = get(groot,'defaultAxesColorOrder');
			options.color = colors(1,:);
			options.scale = 1;
			options = prtUtilAssignStringValuePairs(options,varargin{:});
			
			sz = [1,1]*options.scale;
			l = 2*options.scale;
			holdFlag = ishold;
			h = obj.plotCamera([sz/2,l],obj.R,obj.c,options.color);
			hold on
			obj.plotRect(sz,obj.R,obj.c,'','-','Color',options.color)
			if holdFlag
				hold on
			else
				hold off
			end
			axis equal
		end
		
		%% P
		function val = get.P(obj)
			val = obj.K*obj.R*[eye(3),-obj.c];
    end
		
    %% c
		function obj = set.c(obj,val)
			obj.c = val(:);
		end
		
		%% t
		function val = get.t(obj)
			val = -obj.R*obj.c;
		end
		
		%% f
		function val = get.f(obj)
			val = obj.K(1);
		end
		
		function obj = set.f(obj,val)
			obj.K(1,1) = val;
			obj.K(2,2) = val;
		end
		
		%% x0
		function val = get.x0(obj)
			val = obj.K(1,3);
		end
		
		function obj = set.x0(obj,val)
			obj.K(1,3) = val;
		end
		
		%% y0
		function val = get.y0(obj)
			val = obj.K(2,3);
		end
		
		function obj = set.y0(obj,val)
			obj.K(2,3) = val;
		end
		
	end
	methods (Static)
		function hs = plotCamera(sz,R,center,varargin)
			if nargin==4
				args1 = {'-','Color',varargin{:}};
				args3 = {'--','Color',varargin{:}};
			else
				args1 = varargin;
				args3 = varargin;
			end
			holdFlag = ishold;
			hs(1) = CameraMatrix.plotLine([sz(1),0,0]*R,center,'',args1{:});
			hold on
			hs(2) = CameraMatrix.plotLine([0,sz(2),0]*R,center,'',args1{:});
			hs(3) = CameraMatrix.plotLine([0,0,sz(3)]*R,center,'',args3{:});
			if holdFlag
				hold on
			else
				hold off
			end
		end
		
		function hs = plotLine(v,center,txt,varargin)
			xyz = [center(:)';
				v+center(:)'];
			hs = plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
		end
		
		function plotRect(sz,R,center,txt,varargin)
			xyz = [-sz(2)/2,-sz(1)/2,0;
				sz(2)/2,-sz(1)/2,0;
				sz(2)/2,sz(1)/2,0;
				-sz(2)/2,sz(1)/2,0;
				-sz(2)/2,-sz(1)/2,0]*R;
			X = xyz(1:5,1)+center(1);
			Y = xyz(1:5,2)+center(2);
			Z = xyz(1:5,3)+center(3);
			plot3(X,Y,Z,varargin{:});
			%fill3(X,Y,Z,varargin{end});
			if nargin>3 && ~isempty(txt)
				text(mean(X(:)),mean(Y(:)),mean(Z(:))+1,txt)
			end
		end
		
	end
end