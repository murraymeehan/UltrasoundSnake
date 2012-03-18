X=[    93   118;   120   111;   149   106;   181   107;   219   102;   261   102;   298   102;   325   103;   357   101;   393   106;   431   102;   454   110;   472   145;   475   174;   476   205;   477   256;   477   309;   478   351;   478   383;   476   430;   475   462;   463   493;   435   500;   377   505;   334   513;   280   517;   239   514;   179   514;   134   512;    70   501;    61   472;    58   414;    56   336;    63   277;    67   182;    70   141 ];
I=dicomread('US005.dcm');
figure, imshow(I);
h = impoly(gca, X);
setColor(h,'yellow');
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);