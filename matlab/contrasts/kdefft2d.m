function [xa,ya,z]=kdefft2d(A,inc,h);
A=A';
  r=4;

P=2*inc;


    h1=h;
    h2=h;

[n d]=size(A);    % Establish dimensions of A,
                  % n=number of data points.

x=A(:,1);         % Extract the first column of A into x.
xmin=min(x);      % Find the minimum value of x.
xmax=max(x);      % Find the maximum value of x.
xrange=xmax-xmin; % Find the range of x.

% xa holds the x 'axis' vector, defining a grid of x values where 
% the k.d. function will be evaluated and plotted.

ax=xmin-xrange/2;
bx=xmax+xrange/2;

xa=linspace(ax,bx,inc); 

y=A(:,2);         % Extract the second column of A into y.
ymin=min(y);
ymax=max(y);
yrange=ymax-ymin;

% Establish an 'axis' vector for the y values.

ay=ymin-yrange/2;
by=ymax+yrange/2;

ya=linspace(ay,by,inc);



% Find the binned kernel weights, c.

c=zeros(inc,inc);

deltax=(bx-ax)/(inc-1);
binx=floor((x-ax)/deltax)+1;

deltay=(by-ay)/(inc-1);
biny=floor((y-ay)/deltay)+1;

for i=1:n, % Loop over data points.
  w=(deltax*deltay);

  weight11=(xa(binx(i)+1)-x(i))*(ya(biny(i)+1)-y(i))/w;
  weight12=(x(i)-xa(binx(i)))*(ya(biny(i)+1)-y(i))/w;  
  weight21=(xa(binx(i)+1)-x(i))*(y(i)-ya(biny(i)))/w;
  weight22=(x(i)-xa(binx(i)))*(y(i)-ya(biny(i)))/w;  

  c(binx(i),biny(i))=c(binx(i),biny(i))+weight11;
  c(binx(i)+1,biny(i))=c(binx(i)+1,biny(i))+weight12;
  c(binx(i),biny(i)+1)=c(binx(i),biny(i)+1)+weight21;
  c(binx(i)+1,biny(i)+1)=c(binx(i)+1,biny(i)+1)+weight22;
end;

c=c';

% Obtain the kernel weights.

L1=max(floor(r*h1*(inc-1)/(bx-ax)),floor(r*h2*(inc-1)/(by-ay)));

L=min(L1,inc-1);

kw=zeros(2*inc,2*inc);

[X1,Y1]=meshgrid((bx-ax)*[-L:L]/((inc-1)*h1),(by-ay)*[-L:L]/((inc-1)*h2));

kw((inc-L+1):(inc+L+1),(inc-L+1):(inc+L+1))=(2*pi)^(-2/2)*exp(-0.5*(X1.*X1+Y1.*Y1))/(n*h1*h2);

% Apply 'fftshift' to kw.

kw=fftshift(kw);

% Perform the convolution.

z=real(ifft2(fft2(c,P,P).*fft2(kw)));

z=z(1:inc,1:inc).*(z(1:inc,1:inc)>0);

return

if g~=0,
%  mesh(xa,ya,z);
  surf(xa,ya,z);
  title(['h = [',num2str(h1),',',num2str(h2),']']);   
end;
