class Visualizer{
	protected: bool m[21][3];

	public: virtual int visualize(string dataFileName, string resdir){return 0;}

	protected: double getValue(int type,Cell* cell){
		if(type == 0) return cell->getDensity()/Constants::mp;
		if(type == 1) return cell->getP();
		if(type == 2) return cell->getT();
		if(type == 3) return cell->getV(0)*1.0e-5;
		if(type == 4) return cell->getV(1)*1.0e-5;
		if(type == 5) return cell->getV(2)*1.0e-5;
		if(type == 6) return cell->getV()*1.0e-5;
		if(type >= 7 and type<=9) {
			return ( cell->getSide(type-7, Cell::FORWARD_SIDE)->getB() + cell->getSide(type-7, Cell::BACKWARD_SIDE)->getB() ) / 2.0;
		}
		if(type == 10) {
			double B = 0;
			for(int axis=0;axis<3;axis++){
				B+= pow( (cell->getSide(axis, Cell::FORWARD_SIDE)->getB() 
					+ cell->getSide(axis, Cell::BACKWARD_SIDE)->getB() ) / 2.0, 2);
			}
			return sqrt(B);
		}
		if(type >= 11 and type<=13) {
			double E = 0;
			for(int i=0;i<4;i++){
				E+=cell->getRib(type-11,i)->getE()/4.0;
			}
			return E;
		}
		if(type == 14) {
			double E = 0;
			for(int axis=0;axis<3;axis++){
				E+=pow(  (cell->getRib(axis,0)->getE() 
					+ cell->getRib(axis,1)->getE() 
					+ cell->getRib(axis,2)->getE() 
					+ cell->getRib(axis,3)->getE()
				      )/4.0, 2);
			}
			return sqrt(E);
		}
		if(type == 15) {
			cell->getDivB();
		}
		if(type == 16) {
			return cell->getFullEnergy();
		}
		if(type >= 17 && type <=19) {
			return ( cell->getSide(type-17, Cell::FORWARD_SIDE)->getBbg() + cell->getSide(type-17, Cell::BACKWARD_SIDE)->getBbg() ) / 2.0;
		}
		if(type == 20) {
			double Bbg = 0;
			for(int axis=0;axis<3;axis++){
				Bbg+= pow( (cell->getSide(axis, Cell::FORWARD_SIDE)->getBbg() 
					+ cell->getSide(axis, Cell::BACKWARD_SIDE)->getBbg() ) / 2.0, 2);
			}
			return sqrt(Bbg);
		}
	}

	protected: string getName(int type){
		if(type==0) return "N";
		else if(type==1) return "P";
		else if(type==2) return "T";
		else if(type==3) return "Vx";
		else if(type==4) return "Vy";
		else if(type==5) return "Vz";
		else if(type==6) return "V";
		else if(type==7) return "Bx";
		else if(type==8) return "By";
		else if(type==9) return "Bz";
		else if(type==10) return "B";
		else if(type==11) return "Ex";
		else if(type==12) return "Ey";
		else if(type==13) return "Ez";
		else if(type==14) return "E";
		else if(type==15) return "divB";
		else if(type==16) return "en";
		else if(type==17) return "Bbgx";
		else if(type==18) return "Bbgy";
		else if(type==19) return "Bbgz";
		else if(type==20) return "Bbg";
	}
	protected: int getTypeByName(string name){
		if(name.compare("N") == 0) return 0;
		else if(name.compare("P") == 0) return 1;
		else if(name.compare("T") == 0) return 2;
		else if(name.compare("Vx") == 0) return 3;
		else if(name.compare("Vy") == 0) return 4;
		else if(name.compare("Vz") == 0) return 5;
		else if(name.compare("V") == 0) return 6;
		else if(name.compare("Bx") == 0) return 7;
		else if(name.compare("By") == 0) return 8;
		else if(name.compare("Bz") == 0) return 9;
		else if(name.compare("B") == 0) return 10;
		else if(name.compare("Ex") == 0) return 11;
		else if(name.compare("Ey") == 0) return 12;
		else if(name.compare("Ez") == 0) return 13;
		else if(name.compare("E") == 0) return 14;
		else if(name.compare("divB") == 0) return 15;
		else if(name.compare("en") == 0) return 16;
		else if(name.compare("Bbgx") == 0) return 17;
		else if(name.compare("Bbgy") == 0) return 18;
		else if(name.compare("Bbgz") == 0) return 19;
		else if(name.compare("Bbg") == 0) return 20;
		else return -1;
	}
	protected: string getDimention(int type){
		if(type==0) return "sm^-3";
		else if(type==1) return "din/sm^2";
		else if(type==2) return "K";
		else if(type==3) return "km/s";
		else if(type==4) return "km/s";
		else if(type==5) return "km/s";
		else if(type==6) return "km/s";
		else if(type==7) return "Gs";
		else if(type==8) return "Gs";
		else if(type==9) return "Gs";
		else if(type==10) return "Gs";
		else if(type==11) return "";
		else if(type==12) return "";
		else if(type==13) return "";
		else if(type==14) return "";
		else if(type==15) return "";
		else if(type==16) return "erg";
		else if(type==17) return "Gs";
		else if(type==18) return "Gs";
		else if(type==19) return "Gs";
		else if(type==20) return "Gs";
	}
};
