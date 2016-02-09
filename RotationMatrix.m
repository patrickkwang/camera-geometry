classdef RotationMatrix % < hgsetget
	% class for handling rotation matrices
	
	properties (SetAccess=immutable)
		R
		euler % [phi, theta, psi] = [pan/yaw, tilt/pitch, roll]
		q     % quaternion: q(1)+q(2)i+q(3)j+q(4)k
	end
	
	properties (Dependent)
		phi
		theta
		psi
		
		rodriguez
	end
	
	methods
		function obj = RotationMatrix(string,value)
			switch string
				case 'R'
					obj.R = value;
					obj.euler = obj.matrix2euler(obj.R);
					obj.q = obj.matrix2quaternion(obj.R);
				case 'euler'
					obj.euler = value;
					obj.R = obj.euler2matrix(obj.euler);
					obj.q = obj.matrix2quaternion(obj.R);
				case 'taitbryan'
					obj.R = obj.taitbryan2matrix(value);
					obj.euler = obj.matrix2euler(obj.R);
					obj.q = obj.matrix2quaternion(obj.R);
				case {'quaternion','q'}
					obj.q = value/norm(value); % ensure this is unit length
					obj.R = obj.quaternion2matrix(obj.q);
					obj.euler = obj.matrix2euler(obj.R);
			end
		end
		
		%% phi (pan/yaw)
		function val = get.phi(obj)
			val = obj.euler(1);
		end
		
		%% theta (tilt/pitch)
		function val = get.theta(obj)
			val = obj.euler(2);
		end
		
		%% psi (roll)
		function val = get.psi(obj)
			val = obj.euler(3);
		end
		
		%% Rodriguez representation
		function val = get.rodriguez(obj)
			val = obj.quaternion2rodriguez(obj.q);
		end
		
		%% Bingham parameter
		function A = binghamParameter(obj,x)
			% A = RotationMatrix.binghamParameter(x)
			% produces the parameter A = VZV' of the Bingham distribution
			% whose mode is the rotation quaternion with
			% Z = [x 0 0
			%      0 x 0
			%      0 0 0]
			% (assuming isotropic around mode)
			
			if any(x>0)
				error('x should not contain positive values')
			end
			V = [null(obj.q'),obj.q];
			Z = diag([x.*ones(1,3),0]);
			A = V*Z*V';
			
			% ensure symmetry
			A = (A+A')/2;
		end
		
		%% multiplication
		function val = mtimes(obj1,obj2)
			% called when at least one operand of the (*) operator is a
			% RotationMatrix object. Replaces the object OBJ with OBJ.R.
			
			if isa(obj1,'RotationMatrix')
				R1 = obj1.R;
			else
				R1 = obj1;
			end
			if isa(obj2,'RotationMatrix')
				R2 = obj2.R;
			else
				R2 = obj2;
			end
			val = R1*R2;
		end
		
		%% division
		function val = mrdivide(obj1,obj2)
			% called when at least one operand of the (/) operator is a
			% RotationMatrix object. Replaces the object OBJ with OBJ.R.
			
			if isa(obj1,'RotationMatrix')
				R1 = obj1.R;
			else
				R1 = obj1;
			end
			if isa(obj2,'RotationMatrix')
				R2 = obj2.R;
			else
				R2 = obj2;
			end
			val = R1/R2;
		end
		
		function val = mldivide(obj1,obj2)
			% called when at least one operand of the (\) operator is a
			% RotationMatrix object. Replaces the object OBJ with OBJ.R.
			
			if isa(obj1,'RotationMatrix')
				R1 = obj1.R;
			else
				R1 = obj1;
			end
			if isa(obj2,'RotationMatrix')
				R2 = obj2.R;
			else
				R2 = obj2;
			end
			val = R1\R2;
		end
		
		%% negation
		function val = uminus(obj)
			% called when the operand of the (-) operator is a
			% RotationMatrix object. Negates OBJ.R.
			
			val = RotationMatrix('R',-obj.R);
		end
		
		%% transposition
		function val = ctranspose(obj)
			% called when the operand of the (') operator is a
			% RotationMatrix object. Transposes OBJ.R.
			
			val = RotationMatrix('R',obj.R');
		end
	end
	
	methods (Static)
		%% matrix to Euler angles
		function varargout = matrix2euler(R) % changed this to match euler2matrix
% 			psi = -atan2(R(1,3),R(2,3));
% 			theta = atan2(R(3,3),sqrt(R(3,1)^2+R(3,2)^2))-pi/2;
% 			phi = wrapToPi(atan2(R(3,1),R(3,2))+pi);
			psi = -atan2(R(1,3),R(2,3)); %wrapToPi(atan2(R(1,3),R(2,3))-pi);
			theta = -atan2(sqrt(R(3,1)^2+R(3,2)^2),R(3,3)); %-pi/2;
			phi = atan2(R(3,1),R(3,2))-pi;
			
			if nargin==1
				varargout{1} = [phi, theta, psi];
			elseif nargin==3
				varargout{1} = phi;
				varargout{2} = theta;
				varargout{3} = psi;
			end
		end
		
		%% matrix to quaternion
		function q = matrix2quaternion(R)
			q = 1/2*sqrt(1+R(1,1)+R(2,2)+R(3,3));
      c = 1/(4*max(q(1),eps));
			q = [q;
				c*(R(3,2)-R(2,3));
				c*(R(1,3)-R(3,1));
				c*(R(2,1)-R(1,2))];
      q = q/norm(q);
		end
		
		%% quaternion to matrix
		function R = quaternion2matrix(q)
			R = [q(1).^2+q(2).^2-q(3).^2-q(4).^2,...
				2*q(2).*q(3)-2*q(1).*q(4),...
				2*q(2).*q(4)+2*q(1).*q(3);...
				2*q(2).*q(3)+2*q(1).*q(4),...
				q(1).^2-q(2).^2+q(3).^2-q(4).^2,...
				2*q(3).*q(4)-2*q(1).*q(2);...
				2*q(2).*q(4)-2*q(1).*q(3),...
				2*q(3).*q(4)+2*q(1).*q(2),...
				q(1).^2-q(2).^2-q(3).^2+q(4).^2]...
				./(q(1).^2+q(2).^2+q(3).^2+q(4).^2);
			% wikipedia: Rotation formalisms in three dimensions
			% the following works if q is already normalized
% 			R = [1-2*q(3).^2-2*q(4).^2,...
% 				2*(q(2).*q(3)-q(4).*q(1)),...
% 				2*(q(2).*q(4)+q(3).*q(1));...
% 				2*(q(2).*q(3)+q(4).*q(1)),...
% 				1-2*q(2).^2-2*q(4).^2,...
% 				2*(q(3).*q(4)-q(2).*q(1));...
% 				2*(q(2).*q(4)-q(3).*q(1)),...
% 				2*(q(2).*q(1)+q(3).*q(4)),...
% 				1-2*q(4).^2-2*q(1).^2];
		end
		
		%% quaternion to Rodriguez
		function r = quaternion2rodriguez(q)
			e = q(2:4)/norm(q(2:4));
			theta = 2*acos(q(1));
			r = e*theta;
		end
		
		%% quaternion to Euler angles
		function euler = quaternion2euler(q)
			phi = atan2(q(:,2).*q(:,4)+q(:,3).*q(:,1),-q(:,3).*q(:,4)+q(:,2).*q(:,1));
			theta = acos(-q(:,2).^2-q(:,3).^2+q(:,4).^2+q(:,1).^2);
			psi = atan2(q(:,2).*q(:,4)-q(:,3).*q(:,1),q(:,3).*q(:,4)+q(:,2).*q(:,1));
			euler = [phi,theta,psi];
		end
		
		%% Tait-Bryan angles to matrix
		function R = taitbryan2matrix(taitbryan)
% 			R = RotationMatrix.Rx(taitbryan(3))...
% 				*RotationMatrix.Ry(taitbryan(2))...
% 				*RotationMatrix.Rz(taitbryan(1));
			R = RotationMatrix.Rz(taitbryan(3))...
				*RotationMatrix.Ry(taitbryan(2))...
				*RotationMatrix.Rx(taitbryan(1));
		end
		
		%% Euler angles to matrix
		function R = euler2matrix(euler) % I switched the order!!!
			R = RotationMatrix.Rz(euler(3))...
				*RotationMatrix.Rx(euler(2))...
				*RotationMatrix.Rz(euler(1));
		end
		
		%% Axis-angle to quaternion
		function q = axisangle2quaternion(axisangle)
			axis = axisangle(1:3);
			angle = axisangle(4);
			q = transpose([cos(angle/2),sin(angle/2)*axis]);
		end
		
		%% Rx
		function val = Rx(angle)
			c = cos(angle);
			s = sin(angle);
			val = [1 0 0
				0 c -s
				0 s  c];
		end
		
		%% Ry
		function val = Ry(angle)
			c = cos(angle);
			s = sin(angle);
			val = [c 0 s
				0 1 0
				-s 0 c];
		end
		
		%% Rz
		function val = Rz(angle)
			c = cos(angle);
			s = sin(angle);
			val = [c -s 0
				s  c 0
				0 0 1];
		end
	end
end