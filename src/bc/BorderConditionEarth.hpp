class BorderConditionEarth: public BorderCondition {
	private: TaskData* td;
	private: static const int EARTH_BORDER = 101;
	private: double R;
	private: InitialParametersFunction* ipf;
	private: bool dynamical;

	public: BorderConditionEarth(TaskData* td, InitialParametersFunction* ipf){
		BorderConditionEarth::td = td;
		BorderConditionEarth::ipf = ipf;
		R = td->getLowBorderRadius();
		dynamical = false;
	}

	public: int getMark(){
		return EARTH_BORDER;
	};

	public: void markCell(Cell* cell){
		if(!isBorder(cell)) return;

		cell->setMark(getMark());
		for(int axis=0; axis<3; axis++){
			for(int dir=0; dir<2; dir++){
				Side* side = cell->getSide(axis, dir);
				if(isBorder(side->getBCell()) && isBorder(side->getFCell()))
					side->setMark(getMark());
			}

			for(int irib=0; irib<4; irib++){
				Rib* rib = cell->getRib(axis, irib);
/*
				if((rib->getCell(Cell::ur) == NULL || isBorder(rib->getCell(Cell::ur))) 
					&& (rib->getCell(Cell::ul) == NULL || isBorder(rib->getCell(Cell::ul)))
					&& (rib->getCell(Cell::dl) == NULL || isBorder(rib->getCell(Cell::dl)))
					&& (rib->getCell(Cell::dr) == NULL || isBorder(rib->getCell(Cell::dr)))
				  )
*/
					rib->setMark(getMark());
			}
		}
	}

	public: bool isBorder(Cell* cell){
		int ind = 0;
		double r[3];

		if( ((HierarchyCell*) cell)->isSplitted()) return false;

		for(int axis=0; axis<3; axis++){
			for(int dir=0; dir<2; dir++){
				r[axis] = cell->getR(axis) + (2*dir-1) * cell->getH(axis) / 2.0;

				r[(axis+1)%3] = cell->getR((axis+1)%3) + cell->getH((axis+1)%3) / 2.0;
				r[(axis+2)%3] = cell->getR((axis+2)%3) + cell->getH((axis+2)%3) / 2.0;
				if(pow(r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 0.5)  <= R) ind++;

				r[(axis+1)%3] = cell->getR((axis+1)%3) - cell->getH((axis+1)%3) / 2.0;
				r[(axis+2)%3] = cell->getR((axis+2)%3) + cell->getH((axis+2)%3) / 2.0;
				if(pow(r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 0.5)  <= R) ind++;

				r[(axis+1)%3] = cell->getR((axis+1)%3) + cell->getH((axis+1)%3) / 2.0;
				r[(axis+2)%3] = cell->getR((axis+2)%3) - cell->getH((axis+2)%3) / 2.0;
				if(pow(r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 0.5)  <= R) ind++;

				r[(axis+1)%3] = cell->getR((axis+1)%3) - cell->getH((axis+1)%3) / 2.0;
				r[(axis+2)%3] = cell->getR((axis+2)%3) - cell->getH((axis+2)%3) / 2.0;
				if(pow(r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 0.5)  <= R) ind++;
			}
		}
		if(ind > 0) return true;

		return false;
	};
	public: int getBorderMark(Cell* cell){
		if(!isBorder(cell)) return NULL;
		return getMark();
	};

	public: double getBorderFlowOnSide(Rib* rib, int sideIndex, int flowIndex){
		int vIndex = 0;
		if(sideIndex == Rib::LEFT_SIDE) vIndex = Rib::RIGHT_SIDE;
		else if(sideIndex == Rib::RIGHT_SIDE) vIndex = Rib::LEFT_SIDE;
		else if(sideIndex == Rib::UP_SIDE) vIndex = Rib::DOWN_SIDE;
		else if(sideIndex == Rib::DOWN_SIDE) vIndex = Rib::UP_SIDE;

		Side* vside = rib->getSide(vIndex);
		return vside->getFlow(flowIndex);

//		return 0;
	};

	public: double getBorderU(int axis, int borderDirection, Cell* internalCell, double* inU, int indU){
/*
		double U = inU[indU];
		if(indU == 1+axis) U = - inU[indU];
		return U;
*/

		double rho = ipf->getDensity(internalCell->getCell(axis, borderDirection));
		double P = ipf->getPressure(internalCell->getCell(axis, borderDirection));
//		double P1 = internalCell->getCell(axis, 1-borderDirection)->getP();
//		double P2 = internalCell->getP();
//		double P = 2*P2 - P1;

		if(indU == 0) return rho;
		else if(indU >= 1 && indU <= 3) return 0;
		else if(indU == 4) return P / (Constants::gamma - 1.0) 
					+ Utils::csqr(inU[5], inU[6], inU[7])/8.0/Constants::Pi;
		else if(indU >= 5 && indU <= 7) return inU[indU];
	};

	public: bool isDynamical(){
		return dynamical;
	}
};
