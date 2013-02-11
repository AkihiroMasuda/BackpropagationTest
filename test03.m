## Copyright (C) 2013 akihiro
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## test03

## Author: akihiro <akihiro@MARCH>
## Created: 2013-02-10

# �����o��

function main()
	randn('seed', 1);


	xhani = -10:5:10;
	yhani = -10:5:10;
	[xx yy] = meshgrid(xhani, yhani);
	dat = [xx(:)'; yy(:)'];
	ind1 = find(0.3*(xx(:)).^2>yy(:));
	ind0 = find(0.3*(xx(:)).^2<=yy(:));
	K1 = length(ind1);
	K0 = length(ind0);
	x1 = [xx(ind1)'; yy(ind1)'; ones(1,K1); ones(1,K1); zeros(1,K1)];
	x0 = [xx(ind0)'; yy(ind0)'; ones(1,K0); zeros(1,K0); ones(1,K0)];
	X = [x1 x0];
	
	figure(1); clf;
	hold on;
	plot(x1(1,:), x1(2,:), 'b.');
	plot(x0(1,:), x0(2,:), 'm.');
	title('Answer');
	grid on;
	
	
	
	#keyboard;
	
	
	## �W��������
	%w1 = [1; 1; 1];
	%w2 = [1; 1; 1];
	%h = [1;1;1];
	w1 = rand(3,1)*2-1;
	w2 = rand(3,1)*2-1;
	h1 = rand(3,1)*2-1;
	h2 = rand(3,1)*2-1;
	
	k = 0.3;
	
	N = 1000
	for KK_ = 1:N
		for II_ = 1:size(X,2)
			%% �v�Z
			x = X(:,II_);
			%% �ŏI�o��
			[u g1 g2] = calcans(x, w1, w2, h1, h2, @sigmo);
			%%�w�K�f�[�^�i����l�j
			b = x(end-1:end); 
			%% ���ԁ��o�͑w�̏C����
			e1 = -[g1; g2; 1]*(u(1)-b(1))*u(1)*(1-u(1));
			e2 = -[g1; g2; 1]*(u(2)-b(2))*u(2)*(1-u(2));
			%% ���́����ԑw�̏C���p
			c1 = x(1:3)*(e1(1)*h1(1)+e2(1)*h2(1))*g1*(1-g1);
			c2 = x(1:3)*(e1(2)*h1(2)+e2(2)*h2(2))*g2*(1-g2);
			%% �C��
			w1 = w1 + k*c1;
			w2 = w2 + k*c2;
			h1 = h1 + k*e1;
			h2 = h2 + k*e2;
		end
	
		if mod(KK_,100)==0
			KK_
			%���̃V�X�e���̓��������Ă݂�
%			cu = calcans(X, w1, w2, h1, h2, @sigmo);
			cu = calcans(X, w1, w2, h1, h2, @getclass);
			dis = [X(end-1:end,:); cu]
			err = sum(sum(abs(X(end-1:end,:)- cu),1),2)
		end
	end
	
	ans1 = cu(1) > 0.5;
	ans0 = cu(1) <= 0.5;
	
	% ���ʂ�`��
	xhani2 = -10:1:10; 
	yhani2 = -10:1:10;
	[xx2 yy2] = meshgrid(xhani2, yhani2);
	X2 = [xx2(:)'; yy2(:)'; ones(1,length(xx2(:)))];
	[u] = calcans(X2, w1, w2, h1, h2, @getclass);
	ans1_2 = find(u(1)==1);
	ans0_2 = find(u(1)==0);
	
	figure(2); clf;
	hold on;
	plot(X(1,ans1), X(2,ans1), 'bo');
	plot(X(1,ans0), X(2,ans0), 'mo');
	
	plot(X2(1,ans1_2), X2(2,ans1_2), 'b.');
	plot(X2(1,ans0_2), X2(2,ans0_2), 'm.');
	title('Learning result');
	grid on;
	
	keyboard;
end

function ret = sigmo(d)
	ret = 1./(1+exp(-d));
end

% �j���[�����l�b�g�̓��������߂�
function [cu cg1 cg2] = calcans(X, w1, w2, h1, h2, f)
	cg1 = f(w1'*X(1:3,:));
	cg2 = f(w2'*X(1:3,:));
	cu(1,:) = f(h1'*[cg1; cg2; ones(1,size(X,2))]);
	cu(2,:) = f(h2'*[cg1; cg2; ones(1,size(X,2))]);
end

function ret = getclass(d)
	ind = find(1./(1+exp(d))<0.5);
	ret = zeros(1,length(d));
	ret(ind) = 1;
end
