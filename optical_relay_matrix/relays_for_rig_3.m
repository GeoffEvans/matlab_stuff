%See notes 25/3/14 
% calculation of lens options for top relay telecentric relay
% see notes 24/4/12 gaps between lenses go L11,d1,L12,L21,d2,L22
% fc1 is focal length of compound lens closest to objective, fc2 closest to AOL
% focal lenghts of real lenses go fl11, fl21, fl21, fl22
% compound lens formula 1/f=1/f1+1/f2-d/f1f2
% BFL =f2(d-f1)/(d-(f1+f2)) from Wiki Lens optics

close all
clear
Ltot= (682)*1e-3 %total length m between input and back aperture iris planes
Mag=0.6  % iris to iris magnification
fl11=125e-3;
fl12=2;
fl21=300e-3;
fl22=300e-3;

fc1=Mag/(1+Mag)*Ltot/2
fc2=1/(1+Mag)*Ltot/2
FL12=fc1*fl11/(fl11-fc1);
d1=fl12+fl11-fl11*fl12/fc1;
bfl1=fl12*(d1-fl11)/(d1-(fl11+fl12));
c1=fc1-bfl1;

FL22=fc2*fl21/(fl21-fc2)
d2=fl22+fl21-fl21*fl22/fc2;
bfl2=fl22*(d2-fl21)/(d2-(fl21+fl22));
c2=fc2-bfl2;

iris1_ref=2*fc1;
ref_iris2=2*fc2;
l11=iris1_ref-bfl1-d1;
l22=ref_iris2-bfl2-d2;
Res1=round([l11 d1 bfl1+bfl2 d2 l22]*1000)
Sum1 =sum (Res1)
fs1=[fl11 fl12 fl22 fl21]*1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ltot= 682e-3 
Mag=0.8  % iris to iris magnification
fl11=150e-3;
fl12=2;
fl21=300e-3;
fl22=300e-3;

fc1=Mag/(1+Mag)*Ltot/2
fc2=1/(1+Mag)*Ltot/2
FL12=fc1*fl11/(fl11-fc1)
d1=fl12+fl11-fl11*fl12/fc1;
bfl1=fl12*(d1-fl11)/(d1-(fl11+fl12));
c1=fc1-bfl1;

FL22=fc2*fl21/(fl21-fc2)
d2=fl22+fl21-fl21*fl22/fc2;
bfl2=fl22*(d2-fl21)/(d2-(fl21+fl22));
c2=fc2-bfl2;

iris1_ref=2*fc1;
ref_iris2=2*fc2;
l11=iris1_ref-bfl1-d1;
l22=ref_iris2-bfl2-d2;
Res2=([l11 d1 bfl1+bfl2 d2 l22]*1000)
Sum2 =sum (Res2)
fs2=[fl11 fl12 fl21 fl22]*1000


