%% First Part
clear; 										% below co is characteristic equation coefficients
co = [1 -1.2 .07 .3 -.08];					% characteristic equation in z domain
n=5;
unst=0;										% initially assume the system is stable

for i=1:n
	c(i)=co(n-i+1);							% re-arrangement of coefficients to match the theory
end

s = 0;
for i=1:n
	s = s+c(i);
end

if s <0										% criterion no.2
	unst = 1;
end

if unst == 0								% proceed if stable till now
	s1=0;
	
	for i=1:n
		s1=s1+(-1)^(n-1)*c(i)*(-1)^(i-1);
	end
	
	if s1 < 0								% criterion no.3
		unst = 1;
	end
end

if unst==0									% proceed if stable till now
	i=1;									% starting of table formation
	for j=1:n
		TAB(i,j) = c(j);					% First row
	end
	
	i = i+1;
	for j =1:n
		TAB(i,j)=TAB(i-1,n-j+1); 			% 2nd row is in reverse order of first row
	end

	for i =3:2:2*n-(n-2) 					% Beyond 3rd row for each 2 row computation required only for one row
		nj =n-(i-1)/2; 						% no. of elements in a row
		
		for j =1:nj
			TAB(i,j)=TAB(i-2,1)*TAB(i-1,nj+2-j)-TAB(i-1,1)*TAB(i-2,nj+2-j);
		end
		
		i=i+1; 								% each 2nd row
		for j =1:nj
			TAB(i,j)=TAB(i-1,nj-j+1);
		end
	end										% table completed

	if abs(TAB(1,1))-abs(TAB(2,1))>0 		% first criterion
		unst=1;
	end
	
	for i=3:2:2*n-(n-2)
		if abs(TAB(i,1))<=abs(TAB(i+1,1)) 	% 4th criterion
			unst=1;
		end
	end
end

if unst==1
	fprintf('system is unstable');
else
	fprintf('system is stable');
end

%% Second Part
DT=2e-3;
order =2;
ns = [0 0 1];
ds = [1 2 4];
as = tf(ns,ds);								% s-domain transfer function
zbi=c2d(as,DT,'tustin');					% z-domain transfer function using bilinear
[nzbi,dzbi] = tfdata(zbi,'v');
co = dzbi;									% characteristic equation in z domain
n = order+1;
[ybi]=step(zbi);