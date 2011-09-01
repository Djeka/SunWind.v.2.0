class BorderCondition {
	public: virtual void markCell(Cell* cell){};
	public: virtual bool isBorder(Cell* cell){};
	public: virtual int getBorderMark(Cell* cell){};
	public: virtual int getMark(){};
	public: virtual double getBorderU(int axis, int borderDirection, Cell* internalCell, double* internalU, int indU){};
	public: virtual double getDynamicalBorderU(double time, int axis, int borderDirection, Cell* internalCell, double* internalU, int indU){};
	public: virtual double getBorderFlowOnSide(Rib* rib, int sideIndex, int flowElementIndex){};
	public: virtual bool isDynamical(){};
};
