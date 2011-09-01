#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

class Interpolater{

	public: static void interpolateInsideGDparameters(HierarchyCell* cell){
		double* h = cell->getH();
		double* U0 = getU(cell);
		U0[4] = U0[4] - (U0[1]*U0[1]+U0[2]*U0[2]+U0[3]*U0[3])/2.0/U0[0] - (U0[5]*U0[5]+U0[6]*U0[6]+U0[7]*U0[7])/8.0/Constants::Pi;
		double** gradU = new double*[5];
		for(int i=0; i<5; i++) gradU[i] = new double[3];

		for(int axis = 0; axis < 3; axis++){
			double* Uf = getU(cell->getCell(axis, Cell::FORWARD_SIDE));
			double* Ub = getU(cell->getCell(axis, Cell::BACKWARD_SIDE));

			Uf[4] = Uf[4] - (Uf[1]*Uf[1]+Uf[2]*Uf[2]+Uf[3]*Uf[3])/2.0/Uf[0] - (Uf[5]*Uf[5]+Uf[6]*Uf[6]+Uf[7]*Uf[7])/8.0/Constants::Pi;
			Ub[4] = Ub[4] - (Ub[1]*Ub[1]+Ub[2]*Ub[2]+Ub[3]*Ub[3])/2.0/Ub[0] - (Ub[5]*Ub[5]+Ub[6]*Ub[6]+Ub[7]*Ub[7])/8.0/Constants::Pi;

			double h0 = cell->getH(axis);
			double hf = cell->getCell(axis, Cell::FORWARD_SIDE)->getH(axis);
			double hb = cell->getCell(axis, Cell::BACKWARD_SIDE)->getH(axis);

			for(int i=0; i<5; i++){
				gradU[i][axis] = Utils::minmod((Uf[i] - U0[i])/(h0/2.0 + hf/2.0), (U0[i] - Ub[i])/(h0/2.0 + hf/2.0));
			}

			delete[] Uf;
			delete[] Ub;
		}

		for(int i=0; i<8; i++){
			HierarchyCell* subCell = (HierarchyCell*) cell->getSubCell(i);
			int* p = subCell->getPoint();
			subCell->setDensity(U0[0] 
					+ (2 * p[0] - 1) * gradU[0][0] * h[0] / 4.0
					+ (2 * p[1] - 1) * gradU[0][1] * h[1] / 4.0
					+ (2 * p[2] - 1) * gradU[0][2] * h[2] / 4.0);
			for(int axis=0; axis<3; axis++){
				subCell->setMomentum(U0[1+axis]
					+ (2 * p[0] - 1) * gradU[1+axis][0] * h[0] / 4.0
					+ (2 * p[1] - 1) * gradU[1+axis][1] * h[1] / 4.0
					+ (2 * p[2] - 1) * gradU[1+axis][2] * h[2] / 4.0, axis);
			}
			subCell->setP((Constants::gamma - 1.0) * (U0[4] 
					+ (2 * p[0] - 1) * gradU[4][0] * h[0] / 4.0
					+ (2 * p[1] - 1) * gradU[4][1] * h[1] / 4.0
					+ (2 * p[2] - 1) * gradU[4][2] * h[2] / 4.0));
		}

		delete[] U0;
		for(int i=0; i<5; i++) delete[] gradU[i];
		delete[] gradU;
	};

	private: static double* getU(Cell* cell){
		double* U = new double[8];

		U[0] = cell->getDensity();
		U[1] = cell->getMomentum(0);
		U[2] = cell->getMomentum(1);
		U[3] = cell->getMomentum(2);
		U[4] = cell->getFullEnergy();
		U[5] = (cell->getSide(0, Cell::FORWARD_SIDE)->getB() + cell->getSide(0, Cell::BACKWARD_SIDE)->getB() ) / 2.0;
		U[6] = (cell->getSide(1, Cell::FORWARD_SIDE)->getB() + cell->getSide(1, Cell::BACKWARD_SIDE)->getB() ) / 2.0;
		U[7] = (cell->getSide(2, Cell::FORWARD_SIDE)->getB() + cell->getSide(2, Cell::BACKWARD_SIDE)->getB() ) / 2.0;

		return U;
	}

	public: static void interpolateFieldsOnSubCells(HierarchyCell* cell){
		int x[3];
		Cell* bcell;
		Cell* fcell;

		double** outsideB = interpolateOutsideFields(cell);
		double** indsideB = QPCGAL::calculateQP(outsideB);

		for(int axis=0; axis<3; axis++){
			for(int dir=0; dir<2; dir++){
				x[axis] = dir;
				x[(axis+1)%3]=1; x[(axis+2)%3]=1;
				cell->getSubCell(x)->getSide(axis, dir)->setB(outsideB[axis][4*dir + 0]);

				x[(axis+1)%3]=1; x[(axis+2)%3]=0;
				cell->getSubCell(x)->getSide(axis, dir)->setB(outsideB[axis][4*dir + 1]);

				x[(axis+1)%3]=0; x[(axis+2)%3]=1;
				cell->getSubCell(x)->getSide(axis, dir)->setB(outsideB[axis][4*dir + 2]);

				x[(axis+1)%3]=0; x[(axis+2)%3]=0;
				cell->getSubCell(x)->getSide(axis, dir)->setB(outsideB[axis][4*dir + 3]);
			}

			x[axis] = 1;
			x[(axis+1)%3]=1; x[(axis+2)%3]=1;
			cell->getSubCell(x)->getSide(axis, Cell::BACKWARD_SIDE)->setB(indsideB[axis][0]);
			x[(axis+1)%3]=1; x[(axis+2)%3]=0;
			cell->getSubCell(x)->getSide(axis, Cell::BACKWARD_SIDE)->setB(indsideB[axis][1]);
			x[(axis+1)%3]=0; x[(axis+2)%3]=1;
			cell->getSubCell(x)->getSide(axis, Cell::BACKWARD_SIDE)->setB(indsideB[axis][2]);
			x[(axis+1)%3]=0; x[(axis+2)%3]=0;
			cell->getSubCell(x)->getSide(axis, Cell::BACKWARD_SIDE)->setB(indsideB[axis][3]);
		}

		delete[] indsideB;
		delete[] outsideB;
	}

	private: static double** interpolateOutsideFields(HierarchyCell* cell){
		double** outsideB = new double*[3];
		for(int i=0; i<3; i++) outsideB[i] = new double[8];

		for(int axis=0; axis<3; axis++){
			for(int dir=0; dir<2; dir++){
				double Bff;
				double Bfb;
				double Bbf;
				double Bbb;

				if(cell->getCell(axis, dir) != NULL && ((HierarchyCell*) cell->getCell(axis, dir))->isSplitted()
					&& cell->getDeep() == ((HierarchyCell*) cell->getCell(axis, dir))->getDeep() ) {
					Bff = ((HierarchySide*) cell->getSide(axis, dir))->getSubSide(1, 1)->getB();
					Bfb = ((HierarchySide*) cell->getSide(axis, dir))->getSubSide(1, 0)->getB();
					Bbf = ((HierarchySide*) cell->getSide(axis, dir))->getSubSide(0, 1)->getB();
					Bbb = ((HierarchySide*) cell->getSide(axis, dir))->getSubSide(0, 0)->getB();
				} else {
					HierarchyCell* cell1f = (HierarchyCell*) cell->getCell((axis+1)%3, Cell::FORWARD_SIDE);
					HierarchyCell* cell1b = (HierarchyCell*) cell->getCell((axis+1)%3, Cell::BACKWARD_SIDE);
					HierarchyCell* cell2f = (HierarchyCell*) cell->getCell((axis+2)%3, Cell::FORWARD_SIDE);
					HierarchyCell* cell2b = (HierarchyCell*) cell->getCell((axis+2)%3, Cell::BACKWARD_SIDE);


					double B1f = cell1f->getSide(axis, dir)->getB();
					double B1b = cell1b->getSide(axis, dir)->getB();
					double B2f = cell2f->getSide(axis, dir)->getB();
					double B2b = cell2b->getSide(axis, dir)->getB();

					double h1 = cell->getH((axis+1)%3);
					double h2 = cell->getH((axis+2)%3);
					double B0 = cell->getSide(axis, dir)->getB();

					double Bxx = Utils::minmod(cell1f->getSide(axis, dir)->getB() - cell->getSide(axis, dir)->getB(), 
							 cell->getSide(axis, dir)->getB() - cell1b->getSide(axis, dir)->getB()) / h1;
	//				double Bxx = (B1f - B1b) /2/h1;

					double Byy = Utils::minmod(cell2f->getSide(axis, dir)->getB() - cell->getSide(axis, dir)->getB(), 
							 cell->getSide(axis, dir)->getB() - cell2b->getSide(axis, dir)->getB()) / h2;
	//				double Byy = (B2f - B2b) /2/h2;

					Bff = B0 + Bxx*h1/4 + Byy*h2/4;
					Bfb = B0 + Bxx*h1/4 - Byy*h2/4;
					Bbf = B0 - Bxx*h1/4 + Byy*h2/4;
					Bbb = B0 - Bxx*h1/4 - Byy*h2/4;
				}

				outsideB[axis][4*dir + 0] = Bff;
				outsideB[axis][4*dir + 1] = Bfb;
				outsideB[axis][4*dir + 2] = Bbf;
				outsideB[axis][4*dir + 3] = Bbb;
			}
		}
		return outsideB;
	}	
};
