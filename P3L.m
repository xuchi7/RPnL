function sol = P3L(p1, p2, P1_w, P2_w)

	nl = getProjNorm(p1,p2);
	[Vw Pw] = getVP(P1_w, P2_w);
	
	p1 = [p1; ones(1,3)];
	p2 = [p2; ones(1,3)];
	xs = p1;
	xe = p2;
	
	% cal M
	M1 = kron([1 1 1], [nl nl]');
	M2 = kron([P1_w P2_w]', [1 1 1]);
	A = [M1 .* M2];
	B = [nl nl]';
	BTB = B' * B;
	if rank(BTB) < 3
		sol = [];
		return;
	end
	BA = -inv(BTB)*(B')*A;

	% cal P3L
	nc1= cross(xs(:,1),xe(:,1));  nc1= nc1/norm(nc1);
	nc2= cross(xs(:,2),xe(:,2));  nc2= nc2/norm(nc2);
	nc3= cross(xs(:,3),xe(:,3));  nc3= nc3/norm(nc3);%error:nc3= nc1/norm(nc3)

	Vw1= Vw(:,1);  Vw1= Vw1/norm(Vw1);
	Vw2= Vw(:,2);  Vw2= Vw2/norm(Vw2);
	Vw3= Vw(:,3);  Vw3= Vw3/norm(Vw3);

	Xm = cross(nc1,Vw1); Xm = Xm/norm(Xm); %the X axis of Model frame
	Ym = nc1; %the Y axis of Model frame
	Vm1= cross(Xm,Ym);  Vm1= Vm1/norm(Vm1);%the Z axis of Model frame;

	Rot = [Xm, Ym, Vm1]'; % Rot * [Xm, Ym, Vm1] = I. 

	Vw1new = Rot * Vw1;
	cospsi= [0,0,1]*Vw1new; %the angle between Vm1 and Vw1.
	sinpsi= sqrt(1 - cospsi*cospsi);
	Rmw= [1 0 0; 0 cospsi -sinpsi; 0 sinpsi cospsi];

	nc1 = Rot * nc1; % should be the Y axis, i.e. nc1 = [0,1,0]';
	nc2 = Rot * nc2;
	nc3 = Rot * nc3;

	Vm1 = Rmw*Rot*Vw1; % should be the Z axis, i.e. Vm1 = [0, 0, 1]';
	if 1 - abs(Vm1(3)) > 1e-5
		Rmw = Rmw';
	end

	Vm2 = Rmw*(Rot*Vw2);
	Vm3 = Rmw*(Rot*Vw3);

	A2= Vm2(1); B2= Vm2(2); C2= Vm2(3);
	A3= Vm3(1); B3= Vm3(2); C3= Vm3(3);

	x2= nc2(1); y2= nc2(2); z2= nc2(3);
	x3= nc3(1); y3= nc3(2); z3= nc3(3);

	u11 = -z2*A2*y3*B3 +y2*B2*z3*A3;
	u12 = -y2*A2*z3*B3 +z2*B2*y3*A3;
	u13 = -y2*B2*z3*B3 +z2*B2*y3*B3 +y2*A2*z3*A3 -z2*A2*y3*A3;
	u14 = -y2*B2*x3*C3 +x2*C2*y3*B3;
	u15 = x2*C2*y3*A3 -y2*A2*x3*C3;
	u21 = -x2*A2*y3*B3 +y2*B2*x3*A3;
	u22 = -y2*A2*x3*B3 +x2*B2*y3*A3;
	u23 = x2*B2*y3*B3 -y2*B2*x3*B3 -x2*A2*y3*A3 +y2*A2*x3*A3;
	u24 = y2*B2*z3*C3 -z2*C2*y3*B3;
	u25 = y2*A2*z3*C3 -z2*C2*y3*A3;
	u31 = -x2*A2*z3*A3 +z2*A2*x3*A3;
	u32 = -x2*B2*z3*B3 +z2*B2*x3*B3;
	u33 = x2*A2*z3*B3 -z2*A2*x3*B3 +x2*B2*z3*A3 -z2*B2*x3*A3;
	u34 = z2*A2*z3*C3 +x2*A2*x3*C3 -z2*C2*z3*A3 -x2*C2*x3*A3;
	u35 = -z2*B2*z3*C3 -x2*B2*x3*C3 +z2*C2*z3*B3 +x2*C2*x3*B3;
	u36 = -x2*C2*z3*C3 +z2*C2*x3*C3;

	a4= u11*u11+u12*u12-u13*u13-2*u11*u12+u21*u21+u22*u22-u23*u23...
		-2*u21*u22-u31*u31-u32*u32+u33*u33+2*u31*u32;
	a3= 2*(u11*u14-u13*u15-u12*u14+u21*u24-u23*u25...
		-u22*u24-u31*u34+u33*u35+u32*u34);
	a2= -2*u12*u12+u13*u13+u14*u14-u15*u15+2*u11*u12-2*u22*u22+u23*u23...
		+u24*u24-u25*u25+2*u21*u22+2*u32*u32-u33*u33...
		-u34*u34+u35*u35-2*u31*u32-2*u31*u36+2*u32*u36;%error:-u34*u34+u35+u35-...
	a1= 2*(u12*u14+u13*u15+u22*u24+u23*u25-u32*u34-u33*u35-u34*u36);
	a0= u12*u12+u15*u15+u22*u22+u25*u25-u32*u32-u35*u35-u36*u36-2*u32*u36;
	b3= 2*(u11*u13-u12*u13+u21*u23-u22*u23-u31*u33+u32*u33);
	b2= 2*(u11*u15-u12*u15+u13*u14+u21*u25-u22*u25+u23*u24-u31*u35+u32*u35-u33*u34);
	b1= 2*(u12*u13+u14*u15+u22*u23+u24*u25-u32*u33-u34*u35-u33*u36);
	b0= 2*(u12*u15+u22*u25-u32*u35-u35*u36);

	d0= a0*a0-b0*b0;
	d1= 2*(a0*a1-b0*b1);
	d2= a1*a1+2*a0*a2+b0*b0-b1*b1-2*b0*b2;
	d3= 2*(a0*a3+a1*a2+b0*b1-b1*b2-b0*b3);
	d4= a2*a2+2*a0*a4+2*a1*a3+b1*b1+2*b0*b2-b2*b2-2*b1*b3;
	d5= 2*(a1*a4+a2*a3+b1*b2+b0*b3-b2*b3);
	d6= a3*a3+2*a2*a4+b2*b2-b3*b3+2*b1*b3;
	d7= 2*(a3*a4+b2*b3);
	d8= a4*a4+b3*b3;

	rs= roots([d8 d7 d6 d5 d4 d3 d2 d1 d0]);
	rs= real(rs(abs(imag(rs))<0.01));

	if isempty(rs)
%		disp('no solution');
		sol = [];
		return;
	end

	for i= 1:length(rs)
		cosphi= rs(i);
		
		sign1 = sign(a4 * cosphi^4+a3 * cosphi^3 + a2 * cosphi^2 + a1 * cosphi + a0);
		sign2 = sign(b3 * cosphi^3 + b2 * cosphi^2 + b1 * cosphi + b0);
		
		sinphi= -sign1*sign2*sqrt(abs(1-cosphi*cosphi));
		
		N1= u11*cosphi^2+u12*sinphi^2+u13*cosphi*sinphi+u14*cosphi+u15*sinphi;
		N2= u21*cosphi^2+u22*sinphi^2+u23*cosphi*sinphi+u24*cosphi+u25*sinphi;
		N3= u31*cosphi^2+u32*sinphi^2+u33*cosphi*sinphi+u34*cosphi+u35*sinphi+u36;
		
		costheta = N1/N3;
		sintheta = N2/N3;
		normlize = 1/sqrt(costheta^2+sintheta^2);
		costheta = costheta*normlize;
		sintheta = sintheta*normlize;
		
		R1= [costheta 0 sintheta; 0 1 0; -sintheta 0 costheta];
		R2= [cosphi -sinphi 0; sinphi cosphi 0; 0 0 1];
		
		sol(i).R = (Rot')*R1*R2*Rmw*Rot;
		sol(i).T = BA * sol(i).R(:);
	end
