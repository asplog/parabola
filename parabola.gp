#!/usr/bin/env gnuplot

################################################
# Parabola reflection simulator for gnuplot
#  designed by Kota Nakano
#  last update 07/04 2014
################################################

#Script begins
reset

#definition
F=-0.25; #Focal point
xi=4*F;
phi=pi*6/180.0;	#Diffusion angle [rad]
l=-1.0; #Distance between loudspeaker and parabola [m]
d=0.2; #Radius of loudspeaker [m]
Q = 'Q'; Label for Q
F = 'F'; Label for F
#solve
#P=(px,py) is the intersection point
px=0!=phi?(-(2*tan(phi)*(d/2.0-l*tan(phi))-xi)+sqrt((2*tan(phi)*(d/2.0-l*tan(phi))-xi)**2-4*(tan(phi)**2)*((d/2.0-l*tan(phi))**2)))/(2*(tan(phi)**2)):d**2/4.0/xi;
py=tan(phi)*(px-l)+d/2.0;

##Slope of tangent
tx=2*py/sqrt((xi)**2+(-2*py)**2);
ty=  xi/sqrt((xi)**2+(-2*py)**2);

##Direction of incident wave
ix=cos(phi);
iy=sin(phi);
il=sqrt((px-l)**2+(py-d/2)**2);

##Slope of normal vector
nx=-ty;
ny= tx;

##Slope angles between normal vector and incidence
psi=acos((ix*nx+iy*ny)/sqrt(ix**2+iy**2)/sqrt(nx**2+ny**2));

##Drawing length for normal vector and tangent
nl= py/sin(phi+psi);
tl= py/cos(phi+psi);

##Direction of reflected wave
rx= cos(phi+2*psi)
ry= sin(phi+2*psi)
rl= py/sin(phi+2*psi);

##Point for gathered wave
qx=px-py/tan(phi+2*psi);
qy=0.0;

##Focal point
fx=F;
fy=0.0;

#visualize
clear;
##set paramters
##Use for equivalent ratio
#set size ratio 1.0
#set xrange [-8:2]
#set yrange [-5:5]
##Otherwize
set size noratio
set xrange [l-2*d:2*d];	
set yrange [-py-2*d:py+2*d];
##Labels
set xlabel 'space [m]'
set ylabel 'space [m]'

##Drawing environment
##Draw loudspeaker
set style fill transparent solid 0.5 noborder
set object 1 polygon from \
l, d/2 to \
l, -d/2 to \
l-d/4, -d/4 to \
l-d/2, -d/4 to \
l-d/2, d/4 to \
l-d/4, d/4 to \
l, d/2
set object 1 fc rgb '#000000' fs solid lw 0
##Draw emitted area
set object 2 polygon from \
l, d/2 to \
px, py to \
px,-py to \
l,-d/2 to \
l, d/2
set obj 2 fc rgb '#ffff80' fs solid lw 1;
##Labels
if (0 < strlen(Q))\
set obj 4 rect at qx,qy size char strlen(Q), char 1.6;\
set obj 4 fs empty border -1 front;\
set label 4 at qx,qy Q front center;
if (0 < strlen(F))
set obj 5 rect at fx,fy size char strlen(F), char 1.6;\
set obj 5 fs empty border -1 front;\
set label 5 at fx,fy F front center;
set label 16 at l/2,0 'Emitted area' front centre;

##Draw graphs

set parametric;
set trange [-1:1];

plot py*py*t*t/xi,py*t notitle with filledcurve lt rgb '#ffff80';##Emitted area
replot py*py*t*t/xi,py*t t 'parabola' lt rgb '#ff4040' lw 2;##Parabola
replot px+tx*(tl*1.0*(t-0)),py+ty*(tl*1.0*(t-0)) t 'tangent'   lt rgb '#ff8000' lw 2;##Tangent
replot px+nx*(nl*0.5*(t-1)),py+ny*(nl*0.5*(t-1)) t 'normal'     lt rgb '#c0c000' lw 2;##Normal
replot px+ix*(il*0.5*(t-1)),py+iy*(il*0.5*(t-1)) t 'incidence'  lt rgb '#40ff40' lw 2;##Direction of incident wave
replot px+rx*(rl*0.5*(t-1)),py+ry*(rl*0.5*(t-1)) t 'reflection' lt rgb '#008080' lw 2;##Direction of reflected wave
unset parametric

#Print answers
print(sprintf('F(%f m,%f m)',fx,fy));
print(sprintf('Q(%f m,%f m)',qx,qy));

