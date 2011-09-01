#include <iostream>
   using namespace std;
 
class Matrix 
{
private:
	vector<vector<double> > value;
	int nRow;
	int nCol;
public:
	public: Matrix(){
		nRow = 0;
		nCol = 0;
	};

	public: Matrix(int nRow, int nCol){
		Matrix::nRow = nRow;
		Matrix::nCol = nCol;

		value.resize(nRow);

		for(int i=0; i<nRow; i++){
			vector<double> v1(nCol,0);
			value.at(i) = v1;
		}
	}

	public: Matrix(int nRow, int nCol, double defaultValue){
		Matrix::nRow = nRow;
		Matrix::nCol = nCol;

		value.resize(nRow);

		for(int i=0; i<nRow; i++){
			vector<double> v1(nCol,defaultValue);
			value.at(i) = v1;
		}
	}

	public: int getRows(){
		return nRow;
	}

	public: int getCols(){
		return nCol;
	}

        double &operator() (int p1, int p2);
        Matrix operator+ (Matrix);
        Matrix operator- (Matrix);
        Matrix operator- ();
        Matrix operator* (double Num);
        Matrix operator* (Matrix);
        void operator= (Matrix&);
        vector<double> operator* (vector<double>);
};

ostream& operator <<(ostream& Out, Matrix matrix)
{
        Out <<"\n+--------------------------------------------------------------------------------------------------------------+\n";
        for(int i = 0; i < matrix.getRows(); i++)
        {
                Out << '+';
                for(int j = 0; j < matrix.getCols(); j++)
                {
			Out << "\t";
                        Out << matrix(i,j)<<" |";
                }
                Out << "\n";
        }
        Out <<"+--------------------------------------------------------------------------------------------------------------+\n";
        return Out;
};

double& Matrix ::operator()(int p1, int p2)
{
        if (p1 > nRow || p2 > nCol)
        {       cout << "Error\n" << endl;
//                exit(1);
        }
        return value[p1][p2];
};

Matrix Matrix ::operator+ (Matrix matrix)
{
        Matrix TempMatrix = Matrix(nRow,nCol);
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                        TempMatrix(i,j) = (*this)(i,j) + matrix(i,j);
        return TempMatrix;
};
 
Matrix Matrix ::operator- (Matrix matrix)
{
        Matrix TempMatrix = Matrix(nRow,nCol);
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                        TempMatrix(i,j) = (*this)(i,j) - matrix(i,j);
        return TempMatrix;
};

Matrix Matrix ::operator- ()
{
        Matrix TempMatrix = Matrix(nRow,nCol);
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                        TempMatrix(i,j) = -(*this)(i,j);
        return TempMatrix;
};
 
Matrix Matrix ::operator* (double n)
{
        Matrix TempMatrix = Matrix(nRow,nCol);
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                        TempMatrix(i,j) = (*this)(i,j) * n;
        return TempMatrix;
};

Matrix Matrix ::operator* (Matrix matrix)
{
	Matrix TempMatrix = Matrix(nRow, matrix.getCols());
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < TempMatrix.getCols(); j++)
                {
                        TempMatrix(i,j) = 0;
                        for (int k = 0; k < nCol; k++)
                                TempMatrix(i,j) += (*this)(i,k) * matrix(k,j);
                }

        return TempMatrix;
};

void Matrix ::operator= (Matrix& matrix)
{
	nRow = matrix.getRows();
	nCol = matrix.getCols();

	value.resize(nRow);

	for(int i=0; i<nRow; i++){
		vector<double> v1(nCol,0);
		value.at(i) = v1;
	}

        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                {
			(*this)(i,j) = matrix(i,j);
                }
//        return (*this);
};

vector<double> Matrix ::operator* (vector<double> v)
{
	vector<double> tmpV(nCol, 0);
        for (int i = 0; i < nRow; i++)
                for(int j = 0; j < nCol; j++)
                {
			tmpV[i] += (*this)(i, j) * v[j];
                }

        return tmpV;
};
