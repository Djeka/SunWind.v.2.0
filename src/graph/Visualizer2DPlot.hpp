#include <iostream>
#include <fstream>

using namespace std;

class Visualizer2DPlot: public Visualizer{
	private: CellArray* ca;	
	private: TaskData* td;

	private: double xnorm;
	private: double ynorm;

	public: Visualizer2DPlot(CellArray* ca, TaskData* td, double xnorm, double ynorm){
		Visualizer2DPlot::ca = ca;
		Visualizer2DPlot::td = td;
		Visualizer2DPlot::xnorm = xnorm;
		Visualizer2DPlot::ynorm = ynorm;
	}

	public: int visualize(string dataFileName, string resdir){
		loadVisualizeData(dataFileName);

		for(int type=0;type<21;type++){
			for(int axis=0;axis<3;axis++){
				if(m[type][axis]){visualizePlot(type, axis, resdir);}
			}
		}

		return 1;
	}

	private: void loadVisualizeData(string dataFileName){
		ifstream f(dataFileName.c_str());
		while(!f.eof()){
			string str;
			getline(f, str);			

			string name = str.substr(0,str.find("\t"));
			string yz = str.substr(str.find("\t",0)+1, 1);
			string xz = str.substr( str.find("\t", str.find("\t",0)+1)+1, 1);
			string xy = str.substr( str.find("\t", str.find("\t", str.find("\t",0)+1)+1)+1, 1);

			int type = getTypeByName(name);

			if(type>=0){
				if(yz.compare("Y")==0) m[type][0] = true;
				else m[type][0] = false;
				if(xz.compare("Y")==0) m[type][1] = true;
				else m[type][1] = false;
				if(xy.compare("Y")==0) m[type][2] = true;
				else m[type][2] = false;
			}
		}
		f.close();
	}

	public: void visualizePlot(int type, int axis, string resdir){
		int axisX;
		int axisY;

		if(axis == 0){
			axisX = 1;
			axisY = 2;
		} else if(axis == 1){
			axisX = 0;
			axisY = 2;
		} else if(axis == 2){
			axisX = 0;
			axisY = 1;
		}

		ofstream f((	resdir + "/plot/" + 
				Utils::doubleToString(ca->getTime()) 
				+ "_" 
				+ getName(type) 
				+ "_" 
				+ Utils::getSndAxisNames(axis) 
				+ ".plot").c_str());
		cout << (	"\t" + 
				resdir + "/plot/" + 
				Utils::doubleToString(ca->getTime()) 
				+ "_" 
				+ getName(type) 
				+ "_" 
				+ Utils::getSndAxisNames(axis) 
				+ ".plot") << endl;

		f << " load \"gnuplot/head.gnuplot\"" << endl;
		f << " set output \"" <<resdir<<"/"<<getName(type)<<"/"<<Utils::getSndAxisNames(axis)<<"/"<<ca->getTime()<<"_"<<getName(type)<<"_"
			<<Utils::getSndAxisNames(axis) <<".png\"" << endl;
		f << " set title \"t=" << ca->getTime() << "[sec] ,"<<getName(type)<<"[" << getDimention(type) << "]\"" << endl;
		f << " set xlabel \""<<Utils::getAxisName(axisX)<<"[Re]\"" << endl;
		f << " set ylabel \""<<Utils::getAxisName(axisY)<<"[Re]\"" << endl;
		f << " splot \\" << endl;

		for(Iterator* cit = ca->getIterator(); cit->hasNext();){
			Cell* cell = cit->next();
			if( (
				(abs(cell->getR(axis) - td->getSlice(axis)) < cell->getH(axis)/2.0) ||
				(abs(cell->getR(axis) - td->getSlice(axis)) == cell->getH(axis)/2.0 && (cell->getR(axis) - td->getSlice(axis)) > 0 )
			    )
				&& !((HierarchyCell*) cell)->isSplitted()
				&& cell->getMark() == Cell::INTERNAL_CELL_MARK
			  )
				f << "\"-\" notitle ,\\" << endl;
		}
		f << "\"-\" notitle" << endl;

		int color = 0;
		for(Iterator* cit = ca->getIterator(); cit->hasNext();){
			Cell* cell = cit->next();

			if( (
				(abs(cell->getR(axis) - td->getSlice(axis)) < cell->getH(axis)/2.0) ||
				(abs(cell->getR(axis) - td->getSlice(axis)) == cell->getH(axis)/2.0 && (cell->getR(axis) - td->getSlice(axis)) > 0 )
			    )
				&& !((HierarchyCell*) cell)->isSplitted()
				&& cell->getMark() == Cell::INTERNAL_CELL_MARK
				&& cell->getRadius() > 5 * Constants::Re
			    )
			{
				double value = getValue(type, cell);
//				double value = color;
//				double value = cell->getMark();
//if(value != value) cout << "\t\t\t" << ((HierarchyCell*) cell)->getPath() << "\t" << cell->getDensity() << endl;

				if(type >= 7 && type <= 9 && axis != (type-7) ||
						type >= 17 && type <= 19 && axis != (type-17)){
					double firstIndex;
					if(type >= 7 && type <= 9) firstIndex = 7;
					else if(type >= 17 && type <= 19) firstIndex = 17;
					double faxis = type - firstIndex;
					double ratio = 4.0;
					double* r = cell->getR();
					double h1;
					double h2;

					h1 = cell->getH(axisX) / 2.0;
					h2 = cell->getH(axisY) / 2.0;

//					if(faxis == axisX) h1 = h1 / ratio;
//					else if(faxis == axisY) h2 = h2 / ratio;

					if(faxis == axisX){
						if(firstIndex == 7)	value = cell->getSide(faxis, Cell::BACKWARD_SIDE)->getB();
						else if(firstIndex == 17)	value = cell->getSide(faxis, Cell::BACKWARD_SIDE)->getBbg();
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;
						f << " " << r[axisX] / xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << r[axisX] / xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;

						if(firstIndex == 7)	value = cell->getSide(faxis, Cell::FORWARD_SIDE)->getB();
						else if(firstIndex == 17)	value = cell->getSide(faxis, Cell::FORWARD_SIDE)->getBbg();
						f << " " << r[axisX] / xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << r[axisX] / xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX] + h1) / xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] + h1) / xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;
					} else{
						if(firstIndex == 7)	value = cell->getSide(faxis, Cell::BACKWARD_SIDE)->getB();
						else if(firstIndex == 17)	value = cell->getSide(faxis, Cell::BACKWARD_SIDE)->getBbg();
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY]) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX] + h1)/ xnorm << "\t" << (r[axisY] - h2) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] + h1) / xnorm << "\t" << (r[axisY]) / ynorm << "\t" << value << endl << endl;

						if(firstIndex == 7)	value = cell->getSide(faxis, Cell::FORWARD_SIDE)->getB();
						else if(firstIndex == 17)	value = cell->getSide(faxis, Cell::FORWARD_SIDE)->getBbg();
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] - h1)/ xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX] + h1)/ xnorm << "\t" << (r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] + h1) / xnorm << "\t" << (r[axisY] + h2) / ynorm << "\t" << value << endl << endl;
					}
				} else if(type >= 11 && type <= 13){
					double faxis = type - 11;
					double* r = cell->getR();
					double hx = cell->getH(axisX) / 2.0;
					double hy = cell->getH(axisY) / 2.0;

					if(faxis != axisX && faxis != axisY){
						value = cell->getRib(faxis, Cell::dl)->getE();
						f << " " << (r[axisX] - hx)/ xnorm 	<< "\t" <<(r[axisY]-hy) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] - hx)/ xnorm 	<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX]) / xnorm 		<< "\t" <<(r[axisY]-hy) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl << endl;

						value = cell->getRib(faxis, Cell::ul)->getE();
						f << " " << (r[axisX] - hx)/ xnorm 	<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX] - hx)/ xnorm 	<< "\t" <<(r[axisY]+hy) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX]) / xnorm 		<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]+hy) / ynorm << "\t" << value << endl << endl;

						value = cell->getRib(faxis, Cell::ur)->getE();
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]+hy) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX]+hx) / xnorm 	<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX]+hx)/ xnorm 	<< "\t" <<(r[axisY]+hy) / ynorm << "\t" << value << endl << endl;

						value = cell->getRib(faxis, Cell::dr)->getE();
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]-hy) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX])/ xnorm 		<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl << endl;
						f << " " << (r[axisX]+hx) / xnorm 	<< "\t" <<(r[axisY]-hy) / ynorm << "\t" << value << endl;
						f << " " << (r[axisX]+hx)/ xnorm 	<< "\t" <<(r[axisY]) / ynorm << "\t" << value << endl << endl;
					} else {
						value = getValue(type, cell);
						f << " " << (cell->getR(axisX) - cell->getH(axisX)/2) / xnorm << "\t" << 
							(cell->getR(axisY) - cell->getH(axisY)/2) / ynorm << "\t" <<
							value << endl;

						f << " " << (cell->getR(axisX) - cell->getH(axisX)/2) / xnorm << "\t" << 
							(cell->getR(axisY) + cell->getH(axisY)/2) / ynorm << "\t" <<
							value << endl << endl;

						f << " " << (cell->getR(axisX) + cell->getH(axisX)/2) / xnorm << "\t" << 
							(cell->getR(axisY) - cell->getH(axisY)/2) / ynorm << "\t" <<
							value << endl;

						f << " " << (cell->getR(axisX) + cell->getH(axisX)/2) / xnorm << "\t" << 
							(cell->getR(axisY) + cell->getH(axisY)/2) / ynorm << "\t" <<
							value << endl << endl;
					}
				} else {

					f << " " << (cell->getR(axisX) - cell->getH(axisX)/2) / xnorm << "\t" << 
						(cell->getR(axisY) - cell->getH(axisY)/2) / ynorm << "\t" <<
						value << endl;

					f << " " << (cell->getR(axisX) - cell->getH(axisX)/2) / xnorm << "\t" << 
						(cell->getR(axisY) + cell->getH(axisY)/2) / ynorm << "\t" <<
						value << endl << endl;

					f << " " << (cell->getR(axisX) + cell->getH(axisX)/2) / xnorm << "\t" << 
						(cell->getR(axisY) - cell->getH(axisY)/2) / ynorm << "\t" <<
						value << endl;

					f << " " << (cell->getR(axisX) + cell->getH(axisX)/2) / xnorm << "\t" << 
						(cell->getR(axisY) + cell->getH(axisY)/2) / ynorm << "\t" <<
						value << endl << endl;
				}

				f << " e" << endl;

				if(color == 20) color = 0;
				else color++;
			}
		}
		f.close();
	}
};	
