class Method{
//	public: virtual vector<double> calculateFlow(double Ul[], double Ur[]) {};
	public: virtual vector<double> calculateFlow(double Ul[], double Ur[], double& maxspeed) {};
//	public: virtual vector<double> calculateFlow(double Ul[], double Ur[], double Bdl[], double Bdr[]) {};
//	public: virtual vector<double> calculateFlow(double Ul[], double Ur[], double Bdl[], double Bdr[], double& maxspeed) {};

//	public: virtual vector<double> calculateFlow(double Bx, double Ul[], double Ur[]) {};
	public: virtual vector<double> calculateFlow(double Bx, double Ul[], double Ur[], double& maxspeed) {};
//	public: virtual vector<double> calculateFlow(double Bx, double Ul[], double Ur[], double Bdl[], double Bdr[]) {};
	public: virtual vector<double> calculateFlow(double Bx, double Bdx, double Ul[], double Ur[], double Bdl[], double Bdr[], double& maxspeed) {};

/*
	public: virtual double getMaxSpeed(double Ul[], double Ur[]){return 100;};
	public: virtual double getMaxSpeed(double Ul[], double Ur[], double Bdl[], double Bdr[]){return 100;};

	public: virtual double getMaxSpeed(double Bx, double Ul[], double Ur[]){return 100;};
	public: virtual double getMaxSpeed(double Bx, double Ul[], double Ur[], double Bdl[], double Bdr[]){return 100;};
*/
};
