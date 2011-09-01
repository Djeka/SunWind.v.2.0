class TestRaspadRazriva: public TestObject{
	private: HierarchyCellArray* ca;
	private: TaskData* td;
	private: Method* method;
	private: HierarchyCalculater* calc;
	private: string testName;

	public: TestRaspadRazriva(HierarchyCellArray* ca, TaskData* td, Method* method, HierarchyCalculater* calc){
		TestRaspadRazriva::ca = ca;
		TestRaspadRazriva::td = td;
		TestRaspadRazriva::method = method;
		TestRaspadRazriva::calc = calc;

		testName = "TestRaspadRazriva";
	}

	public: string getTestName(){
		return testName;
	}

	public: int test(){
		clock_t start, end, startTotal;
		double elapsed;
		cout << "[TestRaspadRazriva]\tstart..." << endl;
		cout << "[TestRaspadRazriva]\tInit data is initiating...";
		DataInitiaterTestRaspadRazriva di = DataInitiaterTestRaspadRazriva(ca, td);
		di.initiate();
		cout << "[TestRaspadRazriva]\tdone" << endl;

		CASplitConditionStaticSplit splitCondition = CASplitConditionStaticSplit(ca, td);
		ca->splitCellArray(&splitCondition);

		for(HierarchyCellIterator cit = HierarchyCellIterator(ca); cit.hasNext();){
			HierarchyCell* cell = (HierarchyCell*) cit.next();
			if(cell->getR(0) >= -0.1 && cell->getR(0) <= 0.1 && cell->isSplitted())
				cell->merge();
		}

		cout << "[TestRaspadRazriva]\tVisualizer is initiating...";
		Visualizer1DPlot* vis1D = new Visualizer1DPlot(ca, td, 1, "test/res/1D");
		Visualizer2DPlot* vis2D = new Visualizer2DPlot(ca, td, 1, 1);
		cout << "[TestRaspadRazriva]\tdone" << endl;

		for(int itime=0; itime <= 2000; itime++){
			if(itime%50 == 0)	vis1D->visualize("test/raspad_razriva/visualize1D.data");
//			if(itime%100 == 0)	vis2D->visualize("test/raspad_razriva/visualize2D.data", "test/res");
			ca->setTimeStep(100);

startTotal = clock();

			ca->updateBorderCells();
			calc->calculateFlows();
			ca->showLine(calc, itime);

			ca->updateBorderSides();
			calc->calculateFieldE();
			calc->calculateFieldB();

			ca->calculateIncrements();
			ca->update();

cout.setf(ios::showpoint,ios::floatfield); cout.precision(3);
cout << "Total time: " << ((double)(clock()-startTotal))/CLOCKS_PER_SEC << endl;
cout.setf(ios::scientific,ios::floatfield); cout.precision(10);	
		}

		return 1;
	}
};
