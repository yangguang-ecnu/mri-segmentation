MATITK v.2.4.01 Jan 12 2006
Based on Insight Segmentation and Registration Toolkit (ITK) v.2.4 (http://www.itk.org)


---- readme.txt for windows:
1. copy the DLLs to a desired location
2. set search path of matlab or change current directory to the location of the dll

Example MATLAB commands:
========================
matitk('?')
matitk('f')
matitk('s')
matitk('r')
%------
load mri;
D=squeeze(D);
b=matitk('FCA',[5 0.0625 3],double(D));
c=matitk('SCC',[1.4 10 255],double(b),[],[102 82 25]);subplot(131);
imagesc(squeeze(D(:,:,15)));
axis image;
colormap gray
subplot(132);
imagesc(squeeze(b(:,:,15)));
axis image;
colormap gray
subplot(133);
imagesc(squeeze(c(:,:,15)));
axis image;
colormap gray
%------

Please report any comments/bugs to: Vincent Chu <vwchu@sfu.ca>, Ghassan Hamarneh <hamarneh@cs.sfu.ca>