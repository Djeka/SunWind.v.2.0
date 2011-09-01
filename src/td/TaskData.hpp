#include "iostream"
#include "fstream"
#include <sstream>
#include <cstdlib>

class TaskData{
	private: int n[3];
	private: double h[3];
	private: double Ll[3];
	private: double L[3];
	private: double slice[3];
        private: int NP[3];
        private: double Re;

	private: double splittingDeepMeasure;
	private: double dipole;
	private: int swAxis;
	private: int swDirection;
	private: double swN;
	private: double swT;
	private: double swVx;
	private: double swVy;
	private: double swVz;
	private: double swBx;
	private: double swBy;
	private: double swBz;
	private: double lowBorderRadius;
	private: bool loadDump;
	private: string dumpFile;
	private: string swDataFile;
	private: bool loadSWDataFromFile;

	public: TaskData(string name, double rnorm){
		ifstream f(name.c_str());
		int i = 0;
                Re=rnorm;
		while(!f.eof()){
			string str;
			getline(f, str);
			string svalue = str.substr(0,str.find("\t"));
			double fvalue;
			int ivalue;

			if(str.substr(0,1) == "#") continue;

			if(i == 0){
				dumpFile = svalue;
			} else if(i == 1){
				if(svalue == "yes") loadDump = true;
				else loadDump = false;
			} else if(i>=2 and i<=4) {
				n[i-2] = Utils::getIntFromString(svalue);
			} else if(i>=5 and i<=7){
				Ll[i-5] = Utils::getDoubleFromString(svalue)*rnorm;
			} else if(i>=8 and i<=10){
				L[i-8] = Utils::getDoubleFromString(svalue)*rnorm;
			} else if(i>=11 and i<=13){
				slice[i-11] = Utils::getDoubleFromString(svalue)*rnorm;
			} else if(i==14){
				splittingDeepMeasure = Utils::getDoubleFromString(svalue);
			} else if(i==15){
				dipole = Utils::getDoubleFromString(svalue);
			} else if(i==16){
				lowBorderRadius = Utils::getDoubleFromString(svalue) * Constants::Re;
			} else if(i==17){
				swDataFile = svalue;
			} else if(i==18){
				if(svalue == "yes") loadSWDataFromFile = true;
				else loadSWDataFromFile = false;
			} else if(i==19){
				swAxis = Utils::getIntFromString(svalue);
			} else if(i==20){
				swDirection = Utils::getIntFromString(svalue);
			} else if(i==21){
				 swN = Utils::getDoubleFromString(svalue);
			} else if(i==22){
				 swT = Utils::getDoubleFromString(svalue);
			} else if(i==23){
				 swVx = Utils::getDoubleFromString(svalue) * 1.0e+5;
			} else if(i==24){
				 swVy = Utils::getDoubleFromString(svalue) * 1.0e+5;
			} else if(i==25){
				 swVz = Utils::getDoubleFromString(svalue) * 1.0e+5;
			} else if(i==26){
				 swBx = Utils::getDoubleFromString(svalue);
			} else if(i==27){
				 swBy = Utils::getDoubleFromString(svalue);
			} else if(i==28){
				 swBz = Utils::getDoubleFromString(svalue);
			} else if((i>=29)&&(i<=31)) {
                                sscanf(svalue.c_str(), "%d", &ivalue);
                                NP[i-29]=ivalue;
                        };
			i++;
		}
		f.close();

		for(int i=0; i<3; i++){
			h[i] = L[i] / n[i];
		}

		if(loadDump) loadTaskDataFromDump();
	}

        public: TaskData(TaskData* parenttd,int* X, int* dim) {
              int i;
              for(i=0;i<3;i++) {
                NP[i]=1;
                n[i]=parenttd->getN(i);
                Re=parenttd->getRe();
                L[i]=parenttd->getL(i)/((double)dim[i]);
                Ll[i]=parenttd->getLl(i)+L[i]*((double)X[i]);
                slice[i]=parenttd->getSlice(i);
              };
              splittingDeepMeasure=parenttd->getSplittingDeepMeasure();
	      for(int i=0; i<3; i++){
		h[i] = L[i] / n[i];
	      };
              dipole=parenttd->getEarthMagneticDipole();
        }

	private: void loadTaskDataFromDump(){
		string str;
		ifstream f(this->getDumpFilePath().c_str());

		getline(f, str); // Skip time value
		for(int i=0; i<13; i++){
			getline(f, str);

			if(i>=0 && i<=2) n[i] = Utils::getIntFromString(str);
			else if(i>=3 && i<=5) Ll[i] = Utils::getDoubleFromString(str);
			else if(i>=6 && i<=8) L[i] = Utils::getDoubleFromString(str);
			else if(i>=9 && i<=11) h[i] = Utils::getDoubleFromString(str);
			else if(i==12) lowBorderRadius = Utils::getDoubleFromString(str);
		}
	

		f.close();
	}

	public: int getN(int ind){
		return n[ind];
	}
	public: double getSlice(int ind){
		return slice[ind];
	}
	public: double getH(int ind){
		return h[ind];
	}
	public: double getL(int ind){
		return L[ind];
	}
	public: double getLl(int ind){
		return Ll[ind];
	}
	public: double getSplittingDeepMeasure(){
		return splittingDeepMeasure;
	}
	public: double getEarthMagneticDipole(){
		return dipole;
	}
	public: double getLowBorderRadius(){
		return lowBorderRadius;
	}
	public: int getSwAxis(){
		return swAxis;
	}
	public: int getSwDirection(){
		return swDirection;
	}
	public: double getSwN(){
		return swN;
	}
	public: double getSwT(){
		return swT;
	}
	public: double getSwVx(){
		return swVx;
	}
	public: double getSwVy(){
		return swVy;
	}
	public: double getSwVz(){
		return swVz;
	}
	public: double getSwBx(){
		return swBx;
	}
	public: double getSwBy(){
		return swBy;
	}
	public: double getSwBz(){
		return swBz;
	}
	public: int getNP(int ind){
		return NP[ind];
	}
        public: double getRe(){
                return Re;
        };
	public: bool doLoadDump(){
		return loadDump;
	}
	public: string getDumpFilePath(){
		return dumpFile;
	}
	public: string getSwDataFilePath(){
		return swDataFile;
	}
	public: bool doLoadSWDataFromFile(){
		return loadSWDataFromFile;
	}
};
