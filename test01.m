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

## test01

## Author: akihiro <akihiro@MARCH>
## Created: 2013-02-10

function main()
randn('seed', 1);

K = 10;
x1 = randn(2,K)+[3;3]*ones(1,K);
x2 = randn(2,K)-[3;3]*ones(1,K);
X = [[x1;ones(1,K);zeros(1,K)] [x2;ones(1,K);ones(1,K)]];

## ŒW”‰Šú‰»
w1 = [1; 1; 1];
w2 = [1; 1; 1];
h = [1;1;1];

k = 3;

N = 1000
for KK_ = 1:N
for II_ = 1:size(X,2)
	%% ŒvŽZ
	x = X(:,II_);
	g1 = w1'*x(1:3);
	g1 = sigmo(g1);
	g2 = w2'*x(1:3);
	g2 = sigmo(g2);
	u = h'*[g1;g2;1]; %% ÅIo—Í
	u = sigmo(u);
	b = x(end); %%ŠwKƒf[ƒ^i³‰ð’lj
	%% ’†ŠÔ¨o—Í‘w‚ÌC³—Ê
	e0 =  -1*(u-b)*u*(1-u);
	e1 = -g1*(u-b)*u*(1-u);
	e2 = -g2*(u-b)*u*(1-u);
	%% “ü—Í¨’†ŠÔ‘w‚ÌC³—p
	c1 = x(1:3)*e1*h(1)*g1*(1-g1);
	c2 = x(1:3)*e2*h(2)*g2*(1-g2);
	%% C³
	w1 = w1 + k*c1;
	w2 = w2 + k*c2;
	h = h + k*[e1;e2;e0];
%	keyboard;
end

%¡‚ÌƒVƒXƒeƒ€‚Ì“š‚¦‚ðŒ©‚Ä‚Ý‚é
cg1 = sigmo(w1'*X(1:3,:));
cg2 = sigmo(w2'*X(1:3,:));
cu = sigmo(h'*[cg1; cg2; ones(1,size(X,2))]);

end

keyboard;

figure(1); clf;
hold on;
plot(x1(1,:), x1(2,:), 'b.');
plot(x2(1,:), x2(2,:), 'm.');
grid on;

end

function ret = sigmo(d)
	ret = 1./(1+exp(-d));
end

function ret = getclass(d)
	ind = find(1./(1+exp(d))>0.5);
	ret = zeros(1,length(d));
	ret(ind) = 1;
end

