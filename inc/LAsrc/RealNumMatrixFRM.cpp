#include "../LAinc/RealNumberSystemMatrixFrm.h"

using namespace std;
//using namespace LAP;
using namespace LAP::inSDef;

/*
	some flags:
	//defective		have not finished
	//slow			not efficient
	//unsafe		need improving
*/


namespace LAP
{
	RealSysMatrix createIdentityMatrix(Line order)
	{
		RealSysMatrix tmp_matrix(order, order);
		for (Counting icount = 1; icount <= order; icount++)
		{
			for (Counting jcount = 1; jcount <= order; jcount++)
			{
				if (icount == jcount)
				{
					tmp_matrix(icount, jcount) = 1;
				} 
				else
				{
					tmp_matrix(icount, jcount) = 0;
				}
			}
		}
		return tmp_matrix;
	}
	
	RealSysMatrix createDiagonalMatrix(Line order, RealNumValue lambda)
	{
		RealSysMatrix tmp_matrix(order, order);
		for (Counting icount = 1; icount <= order; icount++)
		{
			for (Counting jcount = 1; jcount <= order; jcount++)
			{
				if (icount == jcount)
				{
					tmp_matrix(icount, jcount) = lambda;
				}
				else
				{
					tmp_matrix(icount, jcount) = 0;
				}
			}
		}
		return tmp_matrix;
	}
	
	RealSysMatrix createCotDiagonalMatrix(Line order, RealNumValue lambda)
	{
		RealSysMatrix tmp_matrix(order, order);
		Line nplusone = 1 + order;
		for (Counting icount = 1; icount <= order; icount++)
		{
			for (Counting jcount = 1; jcount <= order; jcount++)
			{
				if ((icount + jcount) == nplusone)
				{
					tmp_matrix(icount, jcount) = lambda;
				} 
				else
				{
					tmp_matrix(icount, jcount) = 0;
				}
			}
		}
		return tmp_matrix;
	}
	
	RealSysMatrix createCotIdentityMatrix(inSDef::Line order)
	{
		RealSysMatrix tmp_matrix(order, order);
		Line nplusone = 1 + order;
		for (Counting icount = 1; icount <= order; icount++)
		{
			for (Counting jcount = 1; jcount <= order; jcount++)
			{
				if ((icount + jcount) == nplusone)
				{
					tmp_matrix(icount, jcount) = 1;
				}
				else
				{
					tmp_matrix(icount, jcount) = 0;
				}
			}
		}
		return tmp_matrix;
	}

	ostream &operator<<(ostream &stream, RealSysMatrix &matrix)
	{
		for (Counting icount = 1; icount <= matrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= matrix.Column(); jcount++)
			{
				stream << " " << matrix(icount, jcount) << " ";
			}
			stream << endl;
		}
		return stream;
	}

	istream &operator>>(istream &stream, RealSysMatrix &matrix)
	{
		for (Counting icount = 1; icount <= matrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= matrix.Column(); jcount++)
			{
				cout << "input the a" << setw(2) << icount << "," << setw(2) << jcount << "->";
				stream >> matrix(icount, jcount);
			}
		}
		return stream;
	}

	RealSysMatrix operator+(RealSysMatrix & lmatrix, RealSysMatrix & rmatrix)
	{
		//defective
		if (lmatrix.Row() != rmatrix.Row())
		{
			cout << "\nError! different rows..." << endl;
			exit(EXIT_FAILURE);
		}
		if (lmatrix.Column() != rmatrix.Column())
		{
			cout << "\nError! different columns..." << endl;
			exit(EXIT_FAILURE);
		}
		RealSysMatrix result(lmatrix.Row(), lmatrix.Column());
		for (Counting icount = 1; icount <= lmatrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= lmatrix.Column(); jcount++)
			{
				result(icount, jcount) = lmatrix(icount, jcount) + rmatrix(icount, jcount);
			}
		}
		return result;
	}

	RealSysMatrix operator-(RealSysMatrix & lmatrix, RealSysMatrix & rmatrix)
	{
		//defective
		if (lmatrix.Row() != rmatrix.Row())
		{
			cout << "\nError! different rows..." << endl;
			exit(EXIT_FAILURE);
		}
		if (lmatrix.Column() != rmatrix.Column())
		{
			cout << "\nError! different columns..." << endl;
			exit(EXIT_FAILURE);
		}
		RealSysMatrix result(lmatrix.Row(), lmatrix.Column());
		for (Counting icount = 1; icount <= lmatrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= lmatrix.Column(); jcount++)
			{
				result(icount, jcount) = lmatrix(icount, jcount) - rmatrix(icount, jcount);
			}
		}
		return result;
	}

	RealSysMatrix operator*(RealSysMatrix & lmatrix, RealSysMatrix & rmatrix)
	{
		Counting m = lmatrix.Row();
		Counting k = lmatrix.Column();
		Counting n = rmatrix.Column();
		RealNumValue sum = 0;
		RealSysMatrix result(m, n);
		for (Counting i = 1; i <= m; i++)
		{
			for (Counting j = 1; j <= n; j++)
			{
				for (Counting s = 1; s <= k; s++)
				{
					sum += lmatrix(i, s) * rmatrix(s, j);
				}
				result(i, j) = sum;
				sum = 0;
			}
		}
		return result;
	}

	RealSysMatrix &operator+=(RealSysMatrix & lmatrix, RealSysMatrix & rmatrix)
	{
		//defective
		if (lmatrix.Row() != rmatrix.Row())
		{
			cout << "\nError! different rows..." << endl;
			exit(EXIT_FAILURE);
		}
		if (lmatrix.Column() != rmatrix.Column())
		{
			cout << "\nError! different columns..." << endl;
			exit(EXIT_FAILURE);
		}
		for (Counting icount = 1; icount <= lmatrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= lmatrix.Column(); jcount++)
			{
				lmatrix(icount, jcount) += rmatrix(icount, jcount);
			}
		}
		return lmatrix;
	}

	RealSysMatrix & operator-=(RealSysMatrix & lmatrix, RealSysMatrix & rmatrix)
	{
		//defective
		if (lmatrix.Row() != rmatrix.Row())
		{
			cout << "\nError! different rows..." << endl;
			exit(EXIT_FAILURE);
		}
		if (lmatrix.Column() != rmatrix.Column())
		{
			cout << "\nError! different columns..." << endl;
			exit(EXIT_FAILURE);
		}
		for (Counting icount = 1; icount <= lmatrix.Row(); icount++)
		{
			for (Counting jcount = 1; jcount <= lmatrix.Column(); jcount++)
			{
				lmatrix(icount, jcount) -= rmatrix(icount, jcount);
			}
		}
		return lmatrix;
	}

	RealSysMatrix &operator*=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix)
	{
		lmatrix = lmatrix * rmatrix;
		return lmatrix;
	}
}