class BorderConditionSunWind: public BorderCondition {
	private: TaskData* td;
	private: static const int SUNWIND_BORDER = 201;
	private: vector< vector<double> > swData;
	private: bool dynamical;

	public: BorderConditionSunWind(TaskData* td){
		BorderConditionSunWind::td = td;
		dynamical = false;
		if(td->doLoadSWDataFromFile()){
			loadDataFromFile();
			dynamical = true;
		}
	}

	private: void loadDataFromFile(){
		string str;
		int year, month, day, hour, minute, sec, msec;
		double Bx, By, Bz, Vx, Vy, Vz, rho, T;

		int start_time_in_sec = readStartTime();

		ifstream f(td->getSwDataFilePath().c_str());
		while(!f.eof()){
			getline(f, str);
			if(str == "") continue;

			//Format: Year Month Day Hour Min Sec Msec Bx[nT] By[nT] Bz[nT] Vx[km/s] Vy[km/s] Vz[km/s] N[cm^(-3)] T[Kelvin]

			year = Utils::getIntFromString(readValue(str));
			month = Utils::getIntFromString(readValue(str));
			day = Utils::getIntFromString(readValue(str));
			hour = Utils::getIntFromString(readValue(str));
			minute = Utils::getIntFromString(readValue(str));
			sec = Utils::getIntFromString(readValue(str));
			msec = Utils::getIntFromString(readValue(str));

			Bx = Utils::getDoubleFromString(readValue(str)) * 1.0e-5;
			By = Utils::getDoubleFromString(readValue(str)) * 1.0e-5;
			Bz = Utils::getDoubleFromString(readValue(str)) * 1.0e-5;
			Vx = Utils::getDoubleFromString(readValue(str)) * 1.0e+5;
			Vy = Utils::getDoubleFromString(readValue(str)) * 1.0e+5;
			Vz = Utils::getDoubleFromString(readValue(str)) * 1.0e+5;
			rho = Utils::getDoubleFromString(readValue(str)) * Constants::mp;
			T = Utils::getDoubleFromString(readValue(str));

			int time_in_sec = 3600*24*day + 3600 * hour + 60 * minute + sec - start_time_in_sec;

			vector<double> p(9);
			p[0] = time_in_sec;
			p[1] = rho;
			p[2] = Vx;
			p[3] = Vy;
			p[4] = Vz;
			p[5] = T;
			p[6] = Bx;
			p[7] = By;
			p[8] = Bz;

			swData.push_back(p);
		}
		f.close();
	}

	private: int readStartTime(){
		string str;
		int year, month, day, hour, minute, sec, msec;
		ifstream f(td->getSwDataFilePath().c_str());

		getline(f, str);
			year = Utils::getIntFromString(readValue(str));
			month = Utils::getIntFromString(readValue(str));
			day = Utils::getIntFromString(readValue(str));
			hour = Utils::getIntFromString(readValue(str));
			minute = Utils::getIntFromString(readValue(str));
			sec = Utils::getIntFromString(readValue(str));
			msec = Utils::getIntFromString(readValue(str));
		f.close();

		return 3600*24*day + 3600 * hour + 60 * minute + sec;
	}

	private: string readValue(string& str){
		string value = str.substr(0,str.find("\t"));
		str = str.substr(str.find("\t") + 1, str.length() - str.find("\t") + 1);
		return value;
	}

	public: int getMark(){
		return SUNWIND_BORDER;
	};

	public: void markCell(Cell* cell){
		if(!isBorder(cell)) return;

		cell->setMark(getMark());
		for(int axis=0; axis<3; axis++){
			for(int dir=0; dir<2; dir++){
				if(axis == td->getSwAxis() && dir == 1 - td->getSwDirection()) continue;
				cell->getSide(axis, dir)->setMark(getMark());
			}

			for(int irib=0; irib<4; irib++){
				Rib* rib = cell->getRib(axis, irib);

				if(axis == td->getSwAxis()){
					rib->setMark(getMark());
				} else if(axis == (td->getSwAxis() + 1)%3){
					if(td->getSwDirection() == Cell::FORWARD_SIDE){
						if(irib == Cell::ur || irib == Cell::dr)
							rib->setMark(getMark());
					} else if(td->getSwDirection() == Cell::BACKWARD_SIDE){
						if(irib == Cell::ul || irib == Cell::dl)
							rib->setMark(getMark());
					}
				} else if(axis == (td->getSwAxis() + 2)%3){
					if(td->getSwDirection() == Cell::FORWARD_SIDE){
						if(irib == Cell::ur || irib == Cell::ul)
							rib->setMark(getMark());
					} else if(td->getSwDirection() == Cell::BACKWARD_SIDE){
						if(irib == Cell::dr || irib == Cell::dl)
							rib->setMark(getMark());
					}
				}
			}
		}
	}

	public: bool isBorder(Cell* cell){
		if(cell->getCell(td->getSwAxis(), td->getSwDirection()) == NULL 
			&& (cell->getMark() == Cell::INTERNAL_CELL_MARK 
				|| cell->getMark() == Cell::FREE_BORDER_CELL_MARK) ){
			return true;
		}

		return false;
	};
	public: int getBorderMark(Cell* cell){
		if(!isBorder(cell)) return NULL;
		return getMark();
	};

	public: double getBorderFlowOnSide(Rib* rib, int sideIndex, int flowIndex){
		Side* side = rib->getSide(sideIndex);
		return side->getFlow(flowIndex);
	};

	public: double getBorderU(int axis, int dir, Cell* cell, double* inU, int ind){
		double rho = td->getSwN() * Constants::mp;

		if(ind == 0) return rho;
		else if(ind == 1) return rho * td->getSwVx();
		else if(ind == 2) return rho * td->getSwVy();
		else if(ind == 3) return rho * td->getSwVz();
		else if(ind == 4) return rho * Constants::Cv * td->getSwT() + rho * Utils::csqr(td->getSwVx(),td->getSwVy(),td->getSwVz()) / 2.0 
								+ Utils::csqr(td->getSwBx(),td->getSwBy(),td->getSwBz()) / 8.0 / Constants::Pi;
		else if(ind == 5) return td->getSwBx();
		else if(ind == 6) return td->getSwBy();
		else if(ind == 7) return td->getSwBz();
	};

	public: double getDynamicalBorderU(double time, int axis, int dir, Cell* cell, double* inU, int ind){
		vector<double> p;
		vector<double> prev_p = swData[0];
		for(int i=1; i<swData.size(); i++){
			p = swData[i];
			if(time < p[0] && time >= prev_p[0]) break;
			prev_p = p;
		}

		double rho = interpTimeValue(time, prev_p[0], p[0], prev_p[1], p[1]);

		if(ind == 0) return rho;
		else if(ind == 1) return rho * interpTimeValue(time, prev_p[0], p[0], prev_p[2], p[2]);
		else if(ind == 2) return rho * interpTimeValue(time, prev_p[0], p[0], prev_p[3], p[3]);
		else if(ind == 3) return rho * interpTimeValue(time, prev_p[0], p[0], prev_p[4], p[4]);
		else if(ind == 4) {
			double Vx = interpTimeValue(time, prev_p[0], p[0], prev_p[2], p[2]);
			double Vy = interpTimeValue(time, prev_p[0], p[0], prev_p[3], p[3]);
			double Vz = interpTimeValue(time, prev_p[0], p[0], prev_p[4], p[4]);
			double T = interpTimeValue(time, prev_p[0], p[0], prev_p[5], p[5]);
			double Bx = interpTimeValue(time, prev_p[0], p[0], prev_p[6], p[6]);
			double By = interpTimeValue(time, prev_p[0], p[0], prev_p[7], p[7]);
			double Bz = interpTimeValue(time, prev_p[0], p[0], prev_p[8], p[8]);

			return rho * Constants::Cv * T + rho * Utils::csqr(Vx,Vy,Vz) / 2.0 + Utils::csqr(Bx,By,Bz) / 8.0 / Constants::Pi;
		}
		else if(ind == 5) return interpTimeValue(time, prev_p[0], p[0], prev_p[6], p[6]);
		else if(ind == 6) return interpTimeValue(time, prev_p[0], p[0], prev_p[7], p[7]);
		else if(ind == 7) return interpTimeValue(time, prev_p[0], p[0], prev_p[8], p[8]);
	};

	private: double interpTimeValue(double time, double start_time, double end_time, double start_value, double end_value){
		return start_value + (time - start_time) * (end_value - start_value) / (end_time - start_time);
	}

	public: bool isDynamical(){
		return dynamical;
	}
};
