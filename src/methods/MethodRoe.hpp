class MethodRoe: public Method{

	private: static const double sqrt2 = 1.41421356237309514547462185873883e+00;
	private: double ce2,by2,bz2,byz,byz2,cet,cet2,u2,v2,w2,V2,signBx,aA2,as2,af2,by_s,bz_s,Nas,Mas,Naf,Maf,N0,M0,col4,col5,col6,line0,line5,line6;
	private: double gammaMinus1,ry,rz,line0el1,line0el2,line2el,line3el,line5el,line6el,Del;
	private: double sqrt4PiRho;

	private: void calucalate_short_vars_0(){
		gammaMinus1 = Constants::gamma - 1.0;
	}

	private: void calucalate_short_vars_1(double rho,double u,double v,double w,double Bx,double by,double bz){
		u2 = pow(u,2);
		v2 = pow(v,2);
		w2 = pow(w,2);
		V2 = u2 + v2 + w2;
		by2 = pow(by,2);
		bz2 = pow(bz,2);
		byz2 = by2 + bz2;
		byz = sqrt(byz2);
		sqrt4PiRho = sqrt(4.0*Constants::Pi*rho);
		if(Bx >= 0) signBx = 1.0;
		else if(Bx < 0) signBx = -1.0;

		if(byz > 0.0){
			by_s = by/byz;
			bz_s = bz/byz;
		} else{
			by_s = 1/sqrt2;
			bz_s = 1/sqrt2;
		}

		line0 = w*by_s-v*bz_s;
		line5 = bz_s*rho/sqrt4PiRho*signBx;
		line6 = by_s*rho/sqrt4PiRho*signBx;

		col4 = -line0;
		col5 = bz_s*sqrt4PiRho/rho*signBx;
		col6 = by_s*sqrt4PiRho/rho*signBx;

		Del = byz2*rho/4.0/Constants::Pi;
	}

	private: void calucalate_short_vars_2(double rho,double u,double v,double w,double Bx,double by,double bz,double aA, double as, double af,
						double qy, double qz,double alpha){
		double eps = rho/4.0/Constants::Pi;
		Nas=-(cet2+af2)*(cet2-as2)/(eps*cet2*byz2);
		Mas=-af*rho*2.0/sqrt4PiRho*(cet2-as2)*signBx/(eps*cet*byz2);
		Naf=+(cet2+as2)*(af2-cet2)/(eps*cet2*byz2);
		Maf=+as*rho*2.0/sqrt4PiRho*(af2-cet2)*signBx/(eps*cet*byz2);
		N0=-1.0;
		M0= 0.0;

		ry = by/4.0/Constants::Pi + Constants::gamma*qy/rho/(Constants::gamma - 1);
		rz = bz/4.0/Constants::Pi + Constants::gamma*qz/rho/(Constants::gamma - 1);

		line0el1 = 2*Constants::Pi*(v*ry+w*rz);
		line0el2 = by*ry+bz*rz;
		line2el = 2*Constants::Pi*ry;
		line3el = 2*Constants::Pi*rz;
		line5el = rho*by/4.0/Constants::Pi;
		line6el = rho*bz/4.0/Constants::Pi;
	}
////////////////////////////////
	public: MethodRoe(){};

	public: virtual vector<double> calculateFlow(double Bx, double Ul[], double Ur[], double& maxspeed) {
		double Bdl[]={0,0,0};
		double Bdr[]={0,0,0};
		return calculateFlow(Bx, 0, Ul, Ur, Bdl, Bdr, maxspeed);
	};

	public: vector<double> calculateFlow(double Bx, double Bdx, double ULx[], double URx[], double Bdl[], double Bdr[], double& maxspeed) {
		calucalate_short_vars_0();
		double Ul[7], Ur[7];

		if(abs(abs(Bdl[0]) - abs(Bdr[0]))/abs(Bdl[0]) < 1.0e-5) Bdx = Bdl[0];
		if(abs(abs(Bdl[1]) - abs(Bdr[1]))/abs(Bdl[1]) < 1.0e-5) Bdl[1] = (1-1.0e-2)*Bdl[1];
		if(abs(abs(Bdl[2]) - abs(Bdr[2]))/abs(Bdl[2]) < 1.0e-5) Bdl[2] = (1-1.0e-2)*Bdl[2];

/*
cout << "ULx: " <<ULx[0]<<","<<ULx[1]<<","<<ULx[2]<<","<<ULx[3]<<","<<ULx[4]<<","<<ULx[5]<<","<<ULx[6]<<","<<ULx[7]<< endl;
cout << "URx: " <<URx[0]<<","<<URx[1]<<","<<URx[2]<<","<<URx[3]<<","<<URx[4]<<","<<URx[5]<<","<<URx[6]<<","<<URx[7]<< endl;
cout << "Bdl: " <<Bdl[0]<<","<<Bdl[1]<<","<<Bdl[2]<<endl;
cout << "Bdr: " <<Bdr[0]<<","<<Bdr[1]<<","<<Bdr[2]<<endl;
cout << "Bx, Bdx: " << Bx << "\t" << Bdx << endl;
*/

////////////////////////////Add Bd field
		ULx[4] = ULx[4] - Utils::csqr(ULx[5],ULx[6],ULx[7]) / 8.0/Constants::Pi;
		URx[4] = URx[4] - Utils::csqr(URx[5],URx[6],URx[7]) / 8.0/Constants::Pi;
		for(int i=0; i<3; i++){
			ULx[5+i] = ULx[5+i] + Bdl[i];
			URx[5+i] = URx[5+i] + Bdr[i];
		}
		ULx[4] = ULx[4] + Utils::csqr(ULx[5],ULx[6],ULx[7]) / 8.0/Constants::Pi;
		URx[4] = URx[4] + Utils::csqr(URx[5],URx[6],URx[7]) / 8.0/Constants::Pi;
		Bx = Bx + Bdx;
//////////////////////////
		for(int i=0; i<5; i++){Ul[i] = ULx[i];	Ur[i] = URx[i];}
		for(int i=0; i<2; i++){Ul[5+i] = ULx[6+i];	Ur[5+i] = URx[6+i];}

		vector<double> Fl = roe_calc_F_by_U(Ul, Bdl, ULx[5]);
		vector<double> Fr = roe_calc_F_by_U(Ur, Bdr, URx[5]);

		double cel = sqrt(Constants::gamma*gammaMinus1 * 
				( ULx[4]- Utils::csqr(ULx[5],ULx[6],ULx[7])/8.0/Constants::Pi - Utils::csqr(ULx[1],ULx[2],ULx[3])/2.0/ULx[0] )/ULx[0]);
		double cer = sqrt(Constants::gamma*gammaMinus1 * 
				( URx[4]- Utils::csqr(URx[5],URx[6],URx[7])/8.0/Constants::Pi - Utils::csqr(URx[1],URx[2],URx[3])/2.0/URx[0] )/URx[0]);

		Matrix OmR;
		Matrix OmL;
		vector<double> ev(7,0);
		roe_calc_OmegaR_OmegaL_eigenvalues(OmR, OmL, ev, cel, cer, Ul, Ur, Bx);

		double BdyDiv4Pi = 0.5*(Bdl[1]+Bdr[1])/4/Constants::Pi;
		double BdzDiv4Pi = 0.5*(Bdl[2]+Bdr[2])/4/Constants::Pi;

		for(int i=0; i<7; i++)
			for(int j=0; j<7; j++)
				OmL(i,j) = ev[i] * OmL(i,j);
		for(int i=0; i<7; i++)
			OmR(4,i) = OmR(4,i) - BdyDiv4Pi*OmR(5,i) - BdzDiv4Pi*OmR(6,i);
		for(int i=0; i<7; i++){
			OmL(i,5) = OmL(i,5) + BdyDiv4Pi*OmL(i,4);
			OmL(i,6) = OmL(i,6) + BdzDiv4Pi*OmL(i,4);
		}

		Matrix A = OmR*OmL;//TODO: уменьшить количество операций при умножении матриц,тк большинство из них разряжены
////////////////////////////Remove B field
		ULx[4] = ULx[4] - Utils::csqr(ULx[5],ULx[6],ULx[7]) / 8.0/Constants::Pi;
		URx[4] = URx[4] - Utils::csqr(URx[5],URx[6],URx[7]) / 8.0/Constants::Pi;
		for(int i=0; i<3; i++){
			ULx[5+i] = ULx[5+i] - Bdl[i];
			URx[5+i] = URx[5+i] - Bdr[i];
		}
		ULx[4] = ULx[4] + Utils::csqr(ULx[5],ULx[6],ULx[7]) / 8.0/Constants::Pi;
		URx[4] = URx[4] + Utils::csqr(URx[5],URx[6],URx[7]) / 8.0/Constants::Pi;
		Bx = Bx - Bdx;
//////////////////////////
		for(int i=0; i<5; i++){Ul[i] = ULx[i];	Ur[i] = URx[i];}
		for(int i=0; i<2; i++){Ul[5+i] = ULx[6+i];	Ur[5+i] = URx[6+i];}
/////////////////////////
/*
		Matrix L = Matrix(7,7);
		vector<double> dF(7,0);
		vector<double> dUt(7,0);
		for(int i=0; i<7; i++){
			L(i,i) = ev[i];
			dF[i] = Fl[i] - Fr[i];
			dUt[i] = Ul[i] - Ur[i];
		}
		Matrix At = RE*(OmR*L*OmL)*invRU;
		vector<double> dFt = At * dUt;
		for(int i=0; i<7; i++) cout << dFt[i] << "\t" << dF[i] << endl;
*/
/////////////////////////

		vector<double> dU(7,0);
		for(int i=0; i<7; i++) dU[i] = Ur[i] - Ul[i];
		vector<double> F_roe = A * dU;

		double Froe[7];
		for(int i=0; i<7; i++) Froe[i] = 0.5 * (Fl[i] + Fr[i] - F_roe[i]);

		vector<double> F(8,0);

		F[0] = Froe[0];
		F[1] = Froe[1];
		F[2] = Froe[2];
		F[3] = Froe[3];
		F[4] = Froe[4];
		F[6] = Froe[5];
		F[7] = Froe[6];

		maxspeed = 0;
		for(int i=0; i<7; i++) if(ev[i] > maxspeed) maxspeed = ev[i];

		return F;
	};

	private: vector<double> roe_calc_F_by_U(double Uv[], double B0[], double Bx){
		vector<double> F(7,0);

		double pi4 = 4.0 * Constants::Pi;
		double pi8 = 8.0 * Constants::Pi;

		double rho = Uv[0];
		double u = Uv[1]/Uv[0];
		double v = Uv[2]/Uv[0];
		double w = Uv[3]/Uv[0];
		double e = Uv[4];
		double By = Uv[5];
		double Bz = Uv[6];
		double B0x = B0[0];
		double B0y = B0[1];
		double B0z = B0[2];
		double B1x = Bx - B0[0];
		double B1y = By - B0[1];
		double B1z = Bz - B0[2];

		double Ek = rho*Utils::csqr(u,v,w)/2.0;
		double Em = Utils::csqr(Bx,By,Bz)/pi8;
		double P = gammaMinus1*(e - Ek - Em);
		double B1B0 = B1x*B0x + B1y*B0y + B1z*B0z;

		F[0] = rho*u;
		F[1] = rho*u*u + P + Em - Bx*Bx/pi4 - Utils::csqr(B0x,B0y,B0z)/pi8 + B0x*B0x/pi4;
		F[2] = rho*u*v - Bx*By/pi4 + B0x*B0y/pi4;
		F[3] = rho*u*w - Bx*Bz/pi4 + B0x*B0z/pi4;
		F[4] = (e+P + Utils::csqr(B1x,B1y,B1z)/pi8 + B1B0/pi4) * u - Bx/pi4 * (u*B1x+v*B1y+w*B1z);
		F[5] = u*By - v*Bx;
		F[6] = u*Bz - w*Bx;

		return F;
	}

	private: void roe_calc_OmegaR_OmegaL_eigenvalues(Matrix& omegaR, Matrix& omegaL, vector<double>& ev, double ce_l, double ce_r,
													double ULeft[], double URight[], double Bx){
		/*----------------------------------*/
		double theta1 = 1.0;
		double theta2 = 1.0;
		double eta1 = 2.0;
		double eta2 = 2.0;
		/*----------------------------------*/
		vector<double> Sl = roe_get_S_by_U(ULeft,Bx);
		vector<double> Sr = roe_get_S_by_U(URight,Bx);

		double Rl = Sl[0];
		double Ul = Sl[1];
		double Vl = Sl[2];
		double Wl = Sl[3];
		double iPl= Sl[4];
		double Yl = Sl[5];
		double Zl = Sl[6];

		double Rr = Sr[0];
		double Ur = Sr[1];
		double Vr = Sr[2];
		double Wr = Sr[3];
		double iPr= Sr[4];
		double Yr = Sr[5];
		double Zr = Sr[6];
		
		double rho = sqrt(ULeft[0]*URight[0]);
		double u = 0.5 * (Rl*ULeft[1]/ULeft[0] + Rr*URight[1]/URight[0]) / (0.5*(Rl+Rr));
		double v = 0.5 * (Rl*ULeft[2]/ULeft[0] + Rr*URight[2]/URight[0]) / (0.5*(Rl+Rr));
		double w = 0.5 * (Rl*ULeft[3]/ULeft[0] + Rr*URight[3]/URight[0]) / (0.5*(Rl+Rr));
		double iH= 0.5 * ( Rl*(iPl/Rl) + Rr*(iPr/Rr) ) / (0.5*(Rl+Rr));
		double by= 0.5 * (Yl  + Yr) / (0.5*(Rl+Rr));
		double bz= 0.5 * (Zl  + Zr) / (0.5*(Rl+Rr));

		calucalate_short_vars_1(rho,u,v,w,Bx,by,bz);

		double Rav = 0.5*(Rl+Rr);
		double Uav = 0.5*(Ul+Ur);
		double Vav = 0.5*(Vl+Vr);
		double Wav = 0.5*(Wl+Wr);
		double iPav= 0.5*(iPl+iPr);
		double Yav = 0.5*(Yl+Yr);
		double Zav = 0.5*(Zl+Zr);

		double dR = Rl - Rr;
		double dU = Ul - Ur;
		double dV = Vl - Vr;
		double dW = Wl - Wr;
		double diP=iPl - iPr;
		double dY = Yl - Yr;
		double dZ = Zl - Zr;

		double q = (2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma * ( Yav*Yav + Zav*Zav + theta1/4.0*dY*dY + theta2/4.0*dZ*dZ
				+ eta1/4.0/Rav*Yav*dY*dR + eta2/4.0/Rav*Zav*dZ*dR);
		double qy = (2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma * (Yav*Rav + (1.0-theta1)/4.0*dY*dR + 
				(1.0-eta1)/4.0/Rav*dR*dR*Yav);
		double qz = (2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma * (Zav*Rav + (1.0-theta2)/4.0*dZ*dR + 
				(1.0-eta2)/4.0/Rav*dR*dR*Zav);

		double delta = q -  (2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma*byz2*rho;
		double delta_y = qy-(2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma*by*rho;
		double delta_z = qz-(2.0-Constants::gamma)/4.0/Constants::Pi/Constants::gamma*bz*rho;

		double alpha = 0.5*Constants::gamma*(q-qy*by-qz*bz);
		double betta = Constants::gamma*(by*delta_y + bz*delta_z);

		ce2 = gammaMinus1*(iH-0.5*V2-Bx*Bx/4.0/Constants::Pi/rho-byz2*rho/4.0/Constants::Pi);
		double maxce2 = max(pow(ce_l,2), pow(ce_r,2));
		double mince2 = min(pow(ce_l,2), pow(ce_r,2));
		if(ce2 > maxce2) ce2 = maxce2;
		else if(ce2 < mince2) ce2 = mince2;

		cet2 = ce2 + alpha;
		cet = sqrt(cet2);

		vector<double> eigens = roe_calc_eigen_values(u,alpha,by,bz,rho,betta,Bx);
		double af = 0.5 * (eigens[0] - eigens[6]);
		double aA = 0.5 * (eigens[1] - eigens[5]);
		double as = 0.5 * (eigens[2] - eigens[4]);
		for(int i=0; i<7; i++) ev[i] = abs(eigens[i]);

		calucalate_short_vars_2(rho,u,v,w,Bx,by,bz,aA,as,af,qy,qz,alpha);

		Matrix omR = roe_calc_omega_R(as,aA,af,u,v,w,by,bz,iH,q,qy,qz,rho,Bx,alpha,Rav);
		omegaR = omR;
		Matrix omL = roe_calc_omega_L(as,aA,af,u,v,w,by,bz,iH,q,qy,qz,rho,Bx,alpha);
		omegaL = omL;
	}

	private: Matrix roe_calc_A(double rho,double u,double v,double w,double q,double qy,double qz,double Bx,double by,double bz,double iH){
		Matrix A = Matrix(7,7);

		A(0,1) = 1;
		A(1,0) = ((V2-bz*qz-by*qy+q)*Constants::gamma-w2-v2-3*u2)/2.0;
		A(1,1) = 3*u-u*Constants::gamma;
		A(1,2) = v-v*Constants::gamma;
		A(1,3) = w-w*Constants::gamma;
		A(1,4) = Constants::gamma-1;
		A(1,5) = qy*Constants::gamma;
		A(1,6) = qz*Constants::gamma;
		A(2,0) = -u*v;
		A(2,1) = v;
		A(2,2) = u;
		A(2,5) = -Bx/Constants::Pi/4.0;
		A(3,0) = -u*w;
		A(3,1) = w;
		A(3,3) = u;
		A(3,6) = -Bx/Constants::Pi/4.0;
		A(4,0) = ((2*rho*u*w2+2*rho*u*v2+2*rho*u2*u+(-2*bz*qz-2*by*qy
				+2*q)*rho*u)*Constants::Pi*Constants::gamma+(-4*rho*u*iH-2*rho*u*w2-2*rho*u*v2-2*rho*u2*u)*Constants::Pi
				+bz*Bx*rho*w+by*Bx*rho*v+Bx*Bx*u)/(rho*Constants::Pi)/4.0;
		A(4,1) = -(4*rho*u2*Constants::Pi*Constants::gamma+(-4*rho*iH-4*rho*u2)*Constants::Pi+Bx*Bx)/(rho*Constants::Pi)/4.0;
		A(4,2) = -(4*u*v*Constants::Pi*Constants::gamma-4*u*v*Constants::Pi+by*Bx)/Constants::Pi/4.0;
		A(4,3) = -(4*u*w*Constants::Pi*Constants::gamma-4*u*w*Constants::Pi+bz*Bx)/Constants::Pi/4.0;
		A(4,4) = u*Constants::gamma;
		A(4,5) = (4*qy*u*Constants::Pi*Constants::gamma-Bx*v)/Constants::Pi/4.0;
		A(4,6) = (4*qz*u*Constants::Pi*Constants::gamma-Bx*w)/Constants::Pi/4.0;
		A(5,0) = (Bx*v-by*rho*u)/rho;
		A(5,1) = by;
		A(5,2) = -Bx/rho;
		A(5,5) = u;
		A(6,0) = (Bx*w-bz*rho*u)/rho;
		A(6,1) = bz;
		A(6,3) = -Bx/rho;
		A(6,6) = u;

		return A;
	}

	private: vector<double> roe_get_S_by_U(double Values[], double Bx){
		vector<double> S(7,0);

		double R = sqrt(Values[0]);
		double U = R * Values[1]/Values[0];
		double V = R * Values[2]/Values[0];
		double W = R * Values[3]/Values[0];
		double iP= R * ( Values[4] + gammaMinus1*(Values[4] - Utils::csqr(Values[1],Values[2],Values[3])/2/Values[0] 
				- Utils::csqr(Bx,Values[5],Values[6])/8/Constants::Pi ) 
				+ Utils::csqr(Bx,Values[5],Values[6])/8/Constants::Pi ) / Values[0];
		double Y = Values[5] / R;
		double Z = Values[6] / R;

		S[0] = R;
		S[1] = U;
		S[2] = V;
		S[3] = W;
		S[4] = iP;
		S[5] = Y;
		S[6] = Z;

		return S;
	}

	private: vector<double> roe_calc_eigen_values(double u,double alpha,double by,double bz,double rho,double betta,double Bx){
		vector<double> eigens(7,0);

		double aA = abs(Bx) / sqrt4PiRho;
		aA2 = pow(aA,2);

		double lp = 0.5*(ce2 + alpha + aA2 + byz2*rho/4.0/Constants::Pi + betta);
		double lQ = (ce2 + alpha)*aA2;

		as2 = lp-sqrt(lp*lp - lQ);
		double as = sqrt(as2);
		af2 = lp+sqrt(lp*lp - lQ);
		double af = sqrt(af2);

		eigens[0] = u+af;
		eigens[1] = u+aA;
		eigens[2] = u+as;
		eigens[3] = u;
		eigens[4] = u-as;
		eigens[5] = u-aA;
		eigens[6] = u-af;

		return eigens;
	}

	private: Matrix roe_calc_omega_R(double as,double aA,double af,double u,double v,double w,double by,double bz,
						double iH,double q,double qy,double qz,double rho,double Bx,double alpha,double Rav){
		Matrix omegaR = Matrix(7,7);

		vector<double> col = roe_omega_R_calc_col_for_0_as_af(1.0,af,af2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,alpha,rho,af,as);
		for(int irow=0; irow<7; irow++) omegaR(irow,0) = col[irow];

		col = roe_omega_R_calc_col_for_aA(1.0,v,w,by,bz,Bx,rho);
		for(int irow=0; irow<7; irow++) omegaR(irow,1) = col[irow];

		col = roe_omega_R_calc_col_for_0_as_af(1.0,as,as2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,alpha,rho,af,as);
		for(int irow=0; irow<7; irow++) omegaR(irow,2) = col[irow];

		col = roe_omega_R_calc_col_for_0_as_af(1.0,0.0,0.0,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,alpha,rho,af,as);
		for(int irow=0; irow<7; irow++) omegaR(irow,3) = col[irow];

		col = roe_omega_R_calc_col_for_0_as_af(-1.0,as,as2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,alpha,rho,af,as);
		for(int irow=0; irow<7; irow++) omegaR(irow,4) = col[irow];

		col = roe_omega_R_calc_col_for_aA(-1.0,v,w,by,bz,Bx,rho);
		for(int irow=0; irow<7; irow++) omegaR(irow,5) = col[irow];

		col = roe_omega_R_calc_col_for_0_as_af(-1.0,af,af2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,alpha,rho,af,as);
		for(int irow=0; irow<7; irow++) omegaR(irow,6) = col[irow];

		Matrix AU = roe_calc_matrix_AU(u,v,w,iH,q,qy,qz,by,bz,Rav);
		Matrix omR = AU * omegaR;
		return omR;
	}

	private: Matrix roe_calc_matrix_AU(double u,double v,double w,double iH,double q,double qy,double qz,double by,double bz,double Rav){
		Matrix AU = Matrix(7,7);

		AU(0,0) = 2.0;
		AU(1,0) = u;
		AU(1,1) = 1.0;
		AU(2,0) = v;
		AU(2,2) = 1.0;
		AU(3,0) = w;
		AU(3,3) = 1.0;
		AU(4,0) = iH/Constants::gamma-q;
		AU(4,1) = u*gammaMinus1/Constants::gamma;
		AU(4,2) = v*gammaMinus1/Constants::gamma;
		AU(4,3) = w*gammaMinus1/Constants::gamma;
		AU(4,4) = 1.0/Constants::gamma;
		AU(4,5) = -qy;
		AU(4,6) = -qz;
		AU(5,0) = by;
		AU(5,5) = 1.0;
		AU(6,0) = bz;
		AU(6,6) = 1.0;

		return AU;
	}

	private: Matrix roe_calc_omega_L(double as,double aA,double af,double u,double v,double w,double by,double bz,double iH,double q,
											double qy,double qz,double rho,double Bx,double alpha){
		Matrix omegaL = Matrix(7,7);

		vector<double> line = roe_omega_L_calc_col_for_0_as_af(1.0,af,af2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,rho,alpha,af,as);
		for(int icol=0; icol<7; icol++)	omegaL(0,icol) = line[icol];

		line = roe_omega_L_calc_col_for_aA(1.0,v,w,by,bz,Bx,rho);
		for(int icol=0; icol<7; icol++)	omegaL(1,icol) = line[icol];

		line = roe_omega_L_calc_col_for_0_as_af(1.0,as,as2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,rho,alpha,af,as);
		for(int icol=0; icol<7; icol++)	omegaL(2,icol) = line[icol];

		line = roe_omega_L_calc_col_for_0_as_af(1.0,0.0,0.0,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,rho,alpha,af,as);
		for(int icol=0; icol<7; icol++)	omegaL(3,icol) = line[icol];

		line = roe_omega_L_calc_col_for_0_as_af(-1.0,as,as2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,rho,alpha,af,as);
		for(int icol=0; icol<7; icol++)	omegaL(4,icol) = line[icol];

		line = roe_omega_L_calc_col_for_aA(-1.0,v,w,by,bz,Bx,rho);
		for(int icol=0; icol<7; icol++)	omegaL(5,icol) = line[icol];

		line = roe_omega_L_calc_col_for_0_as_af(-1.0,af,af2,u,v,w,by,bz,iH,q,qy,qz,Bx,aA,rho,alpha,af,as);
		for(int icol=0; icol<7; icol++)	omegaL(6,icol) = line[icol];

		vector<double> invD = roe_calc_matrix_invD(as,aA,af,by,bz,rho,alpha,Bx);
		for(int i=0; i<7; i++)
			for(int j=0; j<7; j++)
				omegaL(i,j) = invD[i] * omegaL(i,j);
		return omegaL;
	}

	private: vector<double> roe_calc_matrix_invD(double as,double aA,double af,double by,double bz,double rho,double alpha,double Bx){
		vector<double> invD(7,0);

		double daf = roe_calc_matrix_D_cal_d_form_a(af,af2,by,bz,rho,aA,Bx,alpha,af,as);
		double das = roe_calc_matrix_D_cal_d_form_a(as,as2,by,bz,rho,aA,Bx,alpha,af,as);
		double d0 = roe_calc_matrix_D_cal_d_form_a(0.0,0.0,by,bz,rho,aA,Bx,alpha,af,as);

		invD[0] = 1.0 / daf;
		invD[1] = 1.0 / 2.0;
		invD[2] = 1.0 / das;
		invD[3] =-1.0 / d0;
		invD[4] = invD[2];
		invD[5] = invD[1];
		invD[6] = invD[0];

		return invD;
	}

	private: double roe_calc_matrix_D_cal_d_form_a(double a,double a2,double by,double bz,double rho,double aA,double Bx,double alpha,
																double af,double as){

		double M, N, D;
		roe_calc_N_and_M(M,N, a2);
/*
		if(a == as || a == af){
			if(a == af) alpha_x = sqrt((ce2+alpha-as2)/(af2-as2));
			if(a == as) alpha_x = sqrt((af2-ce2-alpha)/(af2-as2));

			d = -2.0/gammaMinus1*(alpha_x*alpha_x*(a*a+ce2+alpha)+(alpha_x*byz+N)*N*rho/8.0/Constants::Pi);
		} else if(a == 0){
			d = -2.0/gammaMinus1*(a*a+(ce2+alpha)+(1+N)*N/2.0*byz2*rho/4.0/Constants::Pi);
		}
*/
//!!!!!!!!!!!!!!!!!!TEST
		D = -2.0/gammaMinus1*(a*a+cet2 + 0.5*N*(1+N)* Del);
//!!!!!!!!!!!!!!!!!!TEST

		return D;
	}

	private: vector<double> roe_omega_L_calc_col_for_0_as_af(double s,double a,double a2,double u,double v,double w,double by,double bz,double iH,
							double q,double qy,double qz,double Bx,double aA,double rho,double alpha,double af,double as){
		vector<double> line(7,0);
		
		double M,N;
		roe_calc_N_and_M(M,N, a2);
/*
!			if(a.eq.as.or.a.eq.af) then
!				b = sqrt(by**2+bz**2)
!				ry_s = ry/b
!				rz_s = rz/b
!				if(a.eq.af) alpha_x = sqrt((ce**2+alpha-as**2)/(af**2-as**2))
!				if(a.eq.as) alpha_x = sqrt((af**2-ce**2-alpha)/(af**2-as**2))
!
!				L = rho*(alpha_x*b+N)/2.0
!
!				line(1) = alpha_x*((ce**2-a**2+s*a*u)/(gamma-1)-(u**2+v**2+w**2)/2.0)
 !    &						- 2.0*s*Constants::Pi*(v*ry_s+w*rz_s)*M + (by*ry_s+bz*rz_s)*L
!				line(2) = alpha_x*(u-s*a/(gamma-1))
!				line(3) = alpha_x*v+2.0*Constants::Pi*s*ry_s*M
!				line(4) = alpha_x*w+2.0*Constants::Pi*s*rz_s*M
!				line(5) = -alpha_x
!				line(6) = -ry_s*L+alpha_x*rho*by/4.0/Constants::Pi
!				line(7) = -rz_s*L+alpha_x*rho*bz/4.0/Constants::Pi
!			elseif(a.eq.0) then
!				L = (1+N)*rho/2.0
!
!				line(1) = (ce**2-a**2+s*a*u)/(gamma-1) - 0.5*(u**2+v**2+w**2)-2*s*Constants::Pi*(v*ry+w*rz)*M+(by*ry+bz*rz)*L
!				line(2) = u-s*a/(gamma-1)
!				line(3) = v+2*Constants::Pi*M*s*ry
!				line(4) = w+2*Constants::Pi*M*s*rz
!				line(5) = -1.0
!				line(6) = -ry*L+rho*by/4.0/Constants::Pi
!				line(7) = -rz*L+rho*bz/4.0/Constants::Pi
!			else
!				write(*,*) '436: a doesnt equal to as,af and 0'
!			endif
!
*/
//!!!!!!!!!!!!!!!!!TEST
		double L = (1+N)*rho/2.0;

		line[0] = (ce2-a2+s*a*u)/gammaMinus1 - 0.5*V2 - s*M * line0el1 + L * line0el2;
		line[1] = u-s*a/gammaMinus1;
		line[2] = v + s*M * line2el;
		line[3] = w + s*M * line3el;
		line[4] = -1.0;
		line[5] = -ry*L + line5el;
		line[6] = -rz*L + line6el;
//!!!!!!!!!!!!!!!!!!TEST

		return line;
	}

	private: vector<double> roe_omega_R_calc_col_for_0_as_af(double s,double a,double a2,double u,double v,double w,double by,double bz,
				double iH,double q,double qy,double qz,double Bx,double aA,double alpha,double rho,double af,double as){
		vector<double> col(7,0);
		double M,N;
		roe_calc_N_and_M(M,N, a2);

/*			if(a.eq.af.or.a.eq.as) then
!				if(a.eq.af) alpha_x = sqrt((ce**2+alpha-as**2)/(af**2-as**2))
!				if(a.eq.as) alpha_x = sqrt((af**2-ce**2-alpha)/(af**2-as**2))
!
!				b = sqrt(by**2+bz**2)
!
!				if(b.ne.0) then
!					by_s = by/b
!					bz_s = bz/b
!				else
!					by_s = 1.0/sqrt(2.0)
!					bz_s = 1.0/sqrt(2.0)
!				endif
!
!				r5 = alpha_x*(-iH + u**2+v**2+w**2 + 2.0*s*a*u)-(v*by_s+w*bz_s)*s*M + 
!     &						gamma/(gamma-1.0)*((2*a**2-q)*alpha_x - (qy*by_s+qz*bz_s)*N)
!
!				col(1) = alpha_x
!				col(2) = alpha_x*(u+2*s*a)
!				col(3) = alpha_x*v-s*by_s*M
!				col(4) = alpha_x*w-s*bz_s*M
!				col(5) = r5
!				col(6) = by_s*N
!				col(7) = bz_s*N
!			elseif(a.eq.0) then
!				absV = sqrt(v**2+u**2+w**2)
!				r5 = -iH + absV**2 + 2*s*a*u - (v*by+w*bz)*s*M + (2*a**2-q-(qy*by+qz*bz)*N)*gamma / (gamma-1)
!
!				col(1) = 1
!				col(2) = u+2*s*a
!				col(3) = v-s*by*M
!				col(4) = w-s*bz*M
!				col(5) = r5
!				col(6) = by*N
!				col(7) = bz*N
!			else
!				write(*,*) '601: a doesnt equal to as,af and 0'
!			endif
*/
//!!!!!!!!!!!!!!!!!!!!!!TEST
		double r5 = -iH + V2 + s*a* 2*u -s*M* (v*by+w*bz) + (2*a2-q-(qy*by+qz*bz)*N)*Constants::gamma / gammaMinus1;

		col[0] = 1.0;
		col[1] = u+2*s*a;
		col[2] = v-s*by*M;
		col[3] = w-s*bz*M;
		col[4] = r5;
		col[5] = by*N;
		col[6] = bz*N;
//!!!!!!!!!!!!!!!!!!!!!!TEST
		return col;
	}

	private: vector<double> roe_omega_L_calc_col_for_aA(double s,double v,double w,double by,double bz,double Bx,double rho){
		vector<double> line(7,0);

		line[0] = line0;
		line[1] = 0.0;
		line[2] = bz_s;
		line[3] = -by_s;
		line[4] = 0.0;
		line[5] = -s*line5;
		line[6] =  s*line6;

		return line;
	}

	private: vector<double> roe_omega_R_calc_col_for_aA(double s,double v,double w,double by,double bz,double Bx,double rho){
		vector<double> col(7,0);

		col[0] = 0;
		col[1] = 0;
		col[2] = bz_s;
		col[3] = -by_s;
		col[4] = col4;
		col[5] = -s*col5;
		col[6] =  s*col6;

		return col;
	}

	private: void roe_calc_N_and_M(double& M, double& N, double a2){
//		double alpha_f = sqrt((cet2-as2)/(af2-as2));
//		double alpha_s = sqrt((af2-cet2)/(af2-as2));
		

/*
!			if(a.eq.as) then
!				N = -alpha_f*(af**2+ce**2)/ce/sqrt(eps)
!				M = -2.0*alpha_f*af*signBx
!			elseif(a.eq.af) then
!				N = alpha_s*(as**2+ce**2)/ce/sqrt(eps)
!				M = 2.0*alpha_s*as*signBx
!			elseif(a.eq.0) then
!				N = -1.0
!				M = 0.0
!			else
!				write(*,*) '629: a doesnt equal to as,af and 0.'
!				write(*,*) 'a=',a,'af=',af,'as=',as
!			endif
*/
//!!!!!!!!!!!!!!!!!!!!!!TEST
		if(a2 == as2){
			N = Nas;
			M = Mas;
		} else if(a2 == af2){
			N = Naf;
			M = Maf;
		} else if(a2 == 0){
			N = N0;
			M = M0;
		}
//!!!!!!!!!!!!!!!!!!!!!!TEST
	}
};
