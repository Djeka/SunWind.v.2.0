class EquestionsRightPartEarthEquilibriumGravitation: public EquestionsRightPartEquilibrium {
	private: InitialParametersFunction* ipf;
	private: Method* method;

	public: EquestionsRightPartEarthEquilibriumGravitation(InitialParametersFunction* ipf, Method* method){
		EquestionsRightPartEarthEquilibriumGravitation::ipf = ipf;
		EquestionsRightPartEarthEquilibriumGravitation::method = method;
	}

	public: vector<double> getIncrements(double dt, Cell* cell){
		vector<double> dU(8,0);

		for(int axis=0;axis<3;axis++){
			Side* fside = cell->getSide(axis, Cell::FORWARD_SIDE);
			Side* bside = cell->getSide(axis, Cell::BACKWARD_SIDE);

			for(int i=0; i<8; i++) dU[i] += dt * (fside->getEqFlow(i) - bside->getEqFlow(i))/cell->getH(axis);
		}

		return dU;
	};

	public: void calculateEquilibriumFlowOnSide(Side* side){
		if(side == NULL || ((HierarchySide*) side)->isSplitted()) return;

		HierarchyCell* bcell = (HierarchyCell*) side->getBCell();
		HierarchyCell* fcell = (HierarchyCell*) side->getFCell();

		if(bcell != NULL && fcell !=NULL
//				&& bcell->getMark() == Cell::INTERNAL_CELL_MARK 
//				&& fcell->getMark() == Cell::INTERNAL_CELL_MARK
				&& side->getMark() == Side::INTERNAL_SIDE_MARK
		  ){

			side->setEqFlow(calculateFlow(side));
		}
	};

	private: vector<double> calculateFlow(Side* side){
		HierarchyCell* lc = (HierarchyCell*) side->getBCell();
		HierarchyCell* rc = (HierarchyCell*) side->getFCell();
		int axis = side->getAxis();
		vector<double> F(8, 0);
		double Ul[8];
		double Ur[8];
		double Bx = side->getB();

		Ul[0] = ipf->getDensity(lc);
		Ul[1] = ipf->getMomentum(axis, lc);
		Ul[2] = ipf->getMomentum((axis+1)%3, lc);
		Ul[3] = ipf->getMomentum((axis+2)%3, lc);
		Ul[5] = (ipf->getB(lc->getSide(axis, Cell::FORWARD_SIDE)) + ipf->getB(lc->getSide(axis, Cell::BACKWARD_SIDE)) ) / 2.0;
		Ul[6] = (ipf->getB(lc->getSide((axis+1)%3,Cell::FORWARD_SIDE)) + ipf->getB(lc->getSide((axis+1)%3,Cell::BACKWARD_SIDE)) ) / 2.0;
		Ul[7] = (ipf->getB(lc->getSide((axis+2)%3,Cell::FORWARD_SIDE)) + ipf->getB(lc->getSide((axis+2)%3,Cell::BACKWARD_SIDE)) ) / 2.0;
		Ul[4] = ipf->getPressure(lc) / (Constants::gamma - 1.0) 
				+ (pow(Ul[0],2) + pow(Ul[1],2) + pow(Ul[2],2)) / 2.0 / Ul[0]
				+ (pow(Ul[5],2) + pow(Ul[6],2) + pow(Ul[7],2)) / 8.0 / Constants::Pi;

		Ur[0] = ipf->getDensity(rc);
		Ur[1] = ipf->getMomentum(axis, rc);
		Ur[2] = ipf->getMomentum((axis+1)%3, rc);
		Ur[3] = ipf->getMomentum((axis+2)%3, rc);
		Ur[5] = (ipf->getB(rc->getSide(axis, Cell::FORWARD_SIDE)) + ipf->getB(rc->getSide(axis, Cell::BACKWARD_SIDE)) ) / 2.0;
		Ur[6] = (ipf->getB(rc->getSide((axis+1)%3,Cell::FORWARD_SIDE)) + ipf->getB(rc->getSide((axis+1)%3,Cell::BACKWARD_SIDE)) ) / 2.0;
		Ur[7] = (ipf->getB(rc->getSide((axis+2)%3,Cell::FORWARD_SIDE)) + ipf->getB(rc->getSide((axis+2)%3,Cell::BACKWARD_SIDE)) ) / 2.0;
		Ur[4] = ipf->getPressure(rc) / (Constants::gamma - 1.0) 
				+ (pow(Ur[0],2) + pow(Ur[1],2) + pow(Ur[2],2)) / 2.0 / Ur[0]
				+ (pow(Ur[5],2) + pow(Ur[6],2) + pow(Ur[7],2)) / 8.0 / Constants::Pi;

		double maxspeed;
		vector<double> inF = method->calculateFlow(Bx, Ul, Ur, maxspeed);

		F[0] 		= inF[0];
		F[1+axis] 	= inF[1];
		F[1+(axis+1)%3] = inF[2];
		F[1+(axis+2)%3] = inF[3];
		F[4] 		= inF[4];
		F[5+axis] 	= inF[5];
		F[5+(axis+1)%3] = inF[6];
		F[5+(axis+2)%3] = inF[7];

		return F;
	}
};
