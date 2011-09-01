#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h> 

using namespace std;

#include "src/cgal/QPCGAL.hpp"
#include "src/cgal/QPCGALF.hpp"
#include "src/cgal/QPCGALout.hpp"

#include "src/math/Matrix.hpp"

#include "src/cell/SideInt.hpp"
#include "src/cell/RibInt.hpp"
#include "src/cell/CellInt.hpp"
#include "src/rp/EquestionsRightPart.hpp"
#include "src/methods/Method.hpp"
#include "src/utils/Utils.hpp"

#include "src/utils/Constants.hpp"
#include "src/td/TaskData.hpp"
#include "src/cell/Side.hpp"
#include "src/cell/Cell.hpp"
#include "src/cell/Rib.hpp"
#include "src/cell/HierarchySide.hpp"
#include "src/cell/HierarchyRib.hpp"
#include "src/cell/HierarchyCell.hpp"
#include "src/calculaters/CalculaterInt.hpp"
#include "src/bc/BorderCondition.hpp"
#include "src/ca/CellArrayInt.hpp"
#include "src/ca/CellArray.hpp"
#include "src/utils/iterator/Iterator.hpp"
#include "src/ca/CASplitCondition.hpp"

#include "src/di/InitialParametersFunction.hpp"
#include "src/di/InitialParametersFunctionEarth.hpp"

#include "src/rp/EquestionsRightPartEquilibrium.hpp"
#include "src/rp/EquestionsRightPartEarthEquilibriumGravitation.hpp"

#include "src/ca/CASplitConditionStaticSplit.hpp"
#include "src/ca/CASplitConditionEarthDipoleInitial.hpp"
#include "src/utils/iterator/HierarchyCellIterator.hpp"
#include "src/utils/iterator/SubCellIterator.hpp"
#include "src/utils/iterator/BorderCellsIterator.hpp"
#include "src/utils/iterator/BindedCellsIterator.hpp"
#include "src/interp/Interpolater.hpp"
#include "src/ca/HierarchyCellArray.hpp"
#include "src/graph/Visualizer.hpp"
#include "src/graph/Visualizer2DPlot.hpp"
#include "src/graph/Visualizer1DPlot.hpp"
#include "src/di/DataInitiater.hpp"
#include "src/di/DataInitiaterForHierarchyCellArray.hpp"
#include "src/di/DataInitiaterTestRaspadRazriva.hpp"
#include "src/methods/MethodLaksaFridrikhsa.hpp"
#include "src/methods/MethodLaksaFridrikhsaDipole.hpp"
#include "src/methods/MethodRoe.hpp"
#include "src/calculaters/Calculater.hpp"
#include "src/calculaters/HierarchyCalculater.hpp"
#include "src/bc/BorderConditionEarth.hpp"
#include "src/bc/BorderConditionSunWind.hpp"

#include "src/test/TestObject.hpp"
#include "src/test/TestRaspadRazriva.hpp"
#include "src/test/TestHierarchyCellArrayLinks.hpp"

void test(){
	string prefix = "[Testing] ";
	string failedText = " has failed";
	string successText = " passed";

	TaskData td = TaskData("test/raspad_razriva/task.data", 1.0);
	HierarchyCellArray ca = HierarchyCellArray(&td);
	MethodLaksaFridrikhsa method = MethodLaksaFridrikhsa();
//	MethodRoe method = MethodRoe();
	HierarchyCalculater calc = HierarchyCalculater(&ca, &method);

	TestRaspadRazriva test = TestRaspadRazriva(&ca, &td, &method, &calc);
	test.test();
}


int main( void )
{

//	test();
//	return 0;

	clock_t start, end, startTotal;
	double elapsed;
	cout.setf(ios::scientific,ios::floatfield);
	cout.precision(10);

	cout << "Load task.data..." << endl;
	TaskData td = TaskData("task.data", Constants::Re);
	cout << "[done]" << endl;

	InitialParametersFunctionEarth initParams = InitialParametersFunctionEarth(&td);
	BorderConditionEarth earthBC =  BorderConditionEarth(&td, &initParams);
	BorderConditionSunWind swBC =  BorderConditionSunWind(&td);

	HierarchyCellArray ca = HierarchyCellArray(&td);
	ca.addBorderCondition(&earthBC);
	ca.addBorderCondition(&swBC);

	if(!td.doLoadDump()){
		cout << "Splitting..." << endl;
		CASplitConditionEarthDipoleInitial splitCondition = CASplitConditionEarthDipoleInitial(&ca, &td, &initParams);
		ca.splitCellArray(&splitCondition);
		cout << "[done]" << endl;

		ca.defineInternalBorders();

		DataInitiaterForHierarchyCellArray di = DataInitiaterForHierarchyCellArray(&ca, &td, &initParams);
		di.initiate();
	}

	Visualizer1DPlot vis1D = Visualizer1DPlot(&ca, &td, Constants::Re, "res/1D");
	Visualizer2DPlot vis2D = Visualizer2DPlot(&ca, &td, Constants::Re, Constants::Re);

//	MethodLaksaFridrikhsa method = MethodLaksaFridrikhsa();
	MethodRoe method = MethodRoe();
	MethodLaksaFridrikhsaDipole bgmethod = MethodLaksaFridrikhsaDipole();
	HierarchyCalculater calc = HierarchyCalculater(&ca, &method, &bgmethod);

	cout << "Calculating equilibrium right part ..." << endl;
	EquestionsRightPartEarthEquilibriumGravitation rp = EquestionsRightPartEarthEquilibriumGravitation(&initParams, &method);
	ca.calculateEquilibriumFlows(&rp);
	cout << "[done]" << endl;

	if(td.doLoadDump()){
		cout << "Loading dump..." << endl;
		ca.loadCellArray();
		cout << "[done]" << endl;
	}

	for(int itime=0; itime <= 1000000; itime++){
		if(itime%2000 == 0) {
			ca.saveCellArray();
//			vis1D.visualize("visualize1D.data");
			vis2D.visualize("visualize.data", "res");
		}

startTotal = clock();

		ca.setTimeStep(10000);
		ca.updateBorderCells();

		calc.calculateFlows();
		ca.showLine(&calc, itime);

		ca.updateBorderSides();
		calc.calculateFieldE();
		calc.calculateFieldB();

		ca.calculateIncrements();
		ca.calculateRightPart(&rp);
		ca.update();

cout.setf(ios::showpoint,ios::floatfield); cout.precision(3);
cout << " Spent time: " << ((double)(clock()-startTotal))/CLOCKS_PER_SEC << endl;
cout.setf(ios::scientific,ios::floatfield); cout.precision(10);	
	}

	vis2D.visualize("visualize.data", "res");
	cout << "Done!" << endl;
};
