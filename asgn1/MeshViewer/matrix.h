#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cstdlib>
#include <vector>
#include <algorithm>

///////////////////////////////////////
// struct MatrixElement
// elements of Matrix class, represent non-zero entry in a matrix
struct MatrixElement
{
	int row, col;	// row and column of the element
	double value;	// value of the element

	// constructor
	MatrixElement(int r, int c, double v)
		: row(r), col(c), value(v) {}

	// compare function for Matrix::SortMatrix() function
	static bool order (MatrixElement e1, MatrixElement e2)
	{
		if (e1.row < e2.row) return true;
		if (e1.row == e2.row) return (e1.col < e2.col);
		return false;
	}
};

// class declaration
typedef std::vector<MatrixElement> MatrixElementList;

// class Matrix definition
class Matrix 
{
private:
	int m, n;	// size of matrix (m = # of rows, n # of columns)
	MatrixElementList elements;	// list of non-zero entries
	int * rowIndex;				// row indice of non-zero entries

	// fields for CG method
	double* diagInv;
	double* r;
	double* r2;
	double* d;
	double* d2;
	double* q;
	double* s;
	double* s2;

public:
	// constructor & destructor
	Matrix(int m, int n) : m(m), n(n) 
	{
		rowIndex = new int[m+1]; 
		diagInv = new double[m];
		r = new double[m];
		r2 = new double[m];
		d = new double[m];
		d2 = new double[m];
		q = new double[m];
		s = new double[m];
		s2 = new double[m];
	}
	~Matrix() 
	{
		delete[] rowIndex; 
		delete[] r;
		delete[] r2;
		delete[] d;
		delete[] d2;
		delete[] q;
		delete[] s;
		delete[] s2;
	}

	/////////////////////////////////////
	// function AddElement
	// add a new entry into the matrix
	void AddElement(int row, int col, double value)
	{
		elements.push_back(MatrixElement(row, col, value));
	}

	/////////////////////////////////////
	// function SortMatrix
	// sort the matrix elements after you add ALL elements into the matrix
	void SortMatrix()
	{
		sort(elements.begin( ), elements.end( ), MatrixElement::order);

		for (int i=0; i<m+1; i++)
			rowIndex[i] = 0;
		for (int i=0; i<(int)elements.size(); i++)
			rowIndex[elements[i].row + 1] = i + 1;

		for (int i=0; i<m; i++)
			diagInv[i] = 0;
		for (int i=0; i<(int)elements.size(); i++)
			if (elements[i].row == elements[i].col)
				diagInv[elements[i].row] = 1.0 / elements[i].value;
	}


	/////////////////////////////////////
	// function Multiply
	// compute A * xIn = xOut
	// the arrays pointed by xIn and xOut have to be pre-allocated 
	// and have enough space
	void Multiply(double* xIn, double* xOut)
	{
		for (int i=0; i<m; i++)
		{
			double sum = 0;
			for (int j=rowIndex[i]; j<rowIndex[i+1]; j++)
				sum += elements[j].value * xIn[elements[j].col];
			xOut[i] = sum;
		}
	}
	/////////////////////////////////////
	// Multiply PreMultiply
	// compute xIn * A = xOut
	// the arrays pointed by xIn and xOut have to be pre-allocated 
	// and have enough space
	void PreMultiply(double* xIn, double* xOut)
	{
		for (int i=0; i<n; i++) xOut[i] = 0;

		for (int i=0; i<m; i++)
		{
			for (int j=rowIndex[i]; j<rowIndex[i+1]; j++)
				xOut[elements[j].col] += elements[j].value * xIn[i];
		}
	}

	/**********************************************/
	/* function: BCG                              */
	/* description: solve Ax = b for unknowns x   */
	/**********************************************/
	void BCG(double* b, double* x, int maxIter, double tolerance)
	{
		double *temp1 = new double[m];
		for(size_t i=0; i<m; i++) {
			temp1[i] = 0.0;
		}
		double *temp2 = new double[m];
		double *r = new double[m];
		double *d = new double[m];
		double alpha = 0.0;
		double beta = 0.0;
		// For initial 
		this->Multiply(x, temp1);
		for()
		for(size_t i=0; i<m; i++) {
			r[i] = b[i] - temp1[i];
			d[i] = r[i];
		}
		this->Multiply(d, temp2);
		double a_no = 0.0, a_deno = 0.0;
		for(size_t i=0; i<m; i++) {
			a_no += r[i]*r[i];
			a_deno += d[i]*temp2[i];
		}
		alpha = a_no/a_deno;
		// Start to iterate
		for(size_t t=1; t<=maxIter; t++) {
			double b_no = 0.0, b_deno = 0.0;
			for(size_t i=0; i<m; i++) {
				// X(t+1) = X(t) + alpha(t)*d(t)
				x[i] = x[i] + alpha*d[i];
				// beta(t+1) = r(t+1)^2/r(t)^2
				b_deno += r[i]*r[i];
				// r(t+1) = r(t) - alpha(t)Ad(t)
				r[i] = r[i] - alpha*temp2[i];
				b_no += r[i]*r[i];
				beta = b_no/b_deno;
				// d(t+1) = r(t+1) + beta(t+1)d(t)
				d[i] = r[i] + beta*d[i];
			}
			// Ad(t+1)
			this->Multiply(d, temp2);
			// Update alpha
			double a_no = 0.0, a_deno = 0.0;
			for(size_t i=0; i<m; i++) {
				a_no += r[i]*r[i];
				a_deno += d[i]*temp2[i];
			}
			// alpha(t+1) = r(t+1)^2/d(t+1)^T A d(t+1)
			alpha = a_no/a_deno;
		}
		delete temp1;
		delete temp2;
		delete r;
		delete d;
 	}

	// friend operators
	friend ostream & operator<< (ostream & out, const Matrix & r) 
	{
		for (int i=0; i<r.m; i++)
		{
			for(int j=r.rowIndex[i]; j<r.rowIndex[i+1]; j++)
				out << r.elements[j].value << " ";
			out << endl;
		}
		return out;
	}
};

#endif __MATRIX_H__
