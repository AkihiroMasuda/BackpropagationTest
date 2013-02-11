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

## test04

## Author: akihiro <akihiro@MARCH>
## Created: 2013-02-10

# 複数出力（出力層３）

function main()
	randn('seed', 2);


	## 教示データの作成
	xhani = -10:3:10;
	yhani = -10:3:10;
	[xx yy] = meshgrid(xhani, yhani);
	dat = [xx(:)'; yy(:)'];
	ind1 = find( (xx(:)<-3) & (yy(:)>-3));
	ind2 = find(xx(:)> 3);
%	ind3 = find((xx(:)>=-3)&(xx(:)<=3));
	ind3 = (setdiff(1:length(xx(:)),[ind1; ind2]))';
	K1 = length(ind1);
	K2 = length(ind2);
	K3 = length(ind3);
	x1 = [xx(ind1)'; yy(ind1)'; ones(1,K1); ones(1,K1); zeros(1,K1); zeros(1,K1)];
	x2 = [xx(ind2)'; yy(ind2)'; ones(1,K2); zeros(1,K2); ones(1,K2); zeros(1,K2)];
	x3 = [xx(ind3)'; yy(ind3)'; ones(1,K3); zeros(1,K3); zeros(1,K3); ones(1,K3)];
	X = [x1 x2 x3]; %データ [x;y;1;答え1;答え2;答え3]
	
	figure(1); clf;
	hold on;
	plot(x1(1,:), x1(2,:), 'b.');
	plot(x2(1,:), x2(2,:), 'c.');
	plot(x3(1,:), x3(2,:), 'm.');
	title('Answer');
	grid on;
	
	
	%return;
	#keyboard;
	
	
	## 係数初期化
	w1 = rand(3,1)*2-1;
	w2 = rand(3,1)*2-1;
	h1 = rand(3,1)*2-1;
	h2 = rand(3,1)*2-1;
	h3 = rand(3,1)*2-1;
	
	k = 0.4;
	
	N = 1000
	for KK_ = 1:N
		for II_ = 1:size(X,2)
			%% 計算
			x = X(:,II_);
			%% 最終出力
			[u g1 g2] = calcans(x, w1, w2, h1, h2, h3, @sigmo);
			%%学習データ（正解値）
			b = x(end-2:end); 
			%% 中間→出力層の修正量
			e1 = -[g1; g2; 1]*(u(1)-b(1))*u(1)*(1-u(1));
			e2 = -[g1; g2; 1]*(u(2)-b(2))*u(2)*(1-u(2));
			e3 = -[g1; g2; 1]*(u(3)-b(3))*u(3)*(1-u(3));
%			e = -[g1; g2; 1]*((u-b)*u*(1-u))';
			%% 入力→中間層の修正用
			c1 = x(1:3)*(e1(1)*h1(1)+e2(1)*h2(1)+e3(1)*h3(1))*g1*(1-g1);
			c2 = x(1:3)*(e1(2)*h1(2)+e2(2)*h2(2)+e3(2)*h3(2))*g2*(1-g2);
			%% 修正
			w1 = w1 + k*c1;
			w2 = w2 + k*c2;
			h1 = h1 + k*e1;
			h2 = h2 + k*e2;
			h3 = h3 + k*e3;
		end
	
		if mod(KK_,100)==0
			KK_
			%今のシステムの答えを見てみる
			cu = calcans(X, w1, w2, h1, h2, h3, @sigmo);
%			cu = calcans(X, w1, w2, h1, h2, h3, @getclass);
			dis = [X(end-2:end,:); cu]
			err = sum(sum(abs(X(end-2:end,:)- cu),1),2)
		end
	end
	
%	ans1 = cu(1) > 0.5;
%	ans0 = cu(1) <= 0.5;
	
	% 結果を描画
	xhani2 = -10:1:10; 
	yhani2 = -10:1:10;
	[xx2 yy2] = meshgrid(xhani2, yhani2);
	X2 = [xx2(:)'; yy2(:)'; ones(1,length(xx2(:)))];
	[u] = calcans(X2, w1, w2, h1, h2, h3, @getclass);
	ans1_2 = find(u(1,:)==1);
	ans2_2 = find(u(2,:)==1);
	ans3_2 = find(u(3,:)==1);
	
	figure(2); clf;
	hold on;
%	plot(X(1,ans1), X(2,ans1), 'bo');
%	plot(X(1,ans0), X(2,ans0), 'mo');
	
	plot(X2(1,ans1_2), X2(2,ans1_2), 'b.');
	plot(X2(1,ans2_2), X2(2,ans2_2), 'c.');
	plot(X2(1,ans3_2), X2(2,ans3_2), 'm.');
	title('Learning result');
	grid on;
	
	keyboard;
end

function ret = sigmo(d)
	ret = 1./(1+exp(-d));
end

% ニューラルネットの答えを求める
function [cu cg1 cg2] = calcans(X, w1, w2, h1, h2, h3, f)
	cg1 = f(w1'*X(1:3,:));
	cg2 = f(w2'*X(1:3,:));
	cu(1,:) = f(h1'*[cg1; cg2; ones(1,size(X,2))]);
	cu(2,:) = f(h2'*[cg1; cg2; ones(1,size(X,2))]);
	cu(3,:) = f(h3'*[cg1; cg2; ones(1,size(X,2))]);
end

function ret = getclass(d)
	ind = find(1./(1+exp(d))<0.5);
	ret = zeros(1,length(d));
	ret(ind) = 1;
end

