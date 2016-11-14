#ifndef _DETERMINANT_GERANAL_H
#define _DETERMINANT_GERANAL_H

#include "DeterminantDeclaration.h"
#include <string>

namespace LAP
{
	
	inline long double calcDiagonalDet(Determinant *diag);

	long double calcVandermondeDet(Determinant *vanDet);
#define CalcUpperTriangularDet(x) calcDiagonalDet(x)
#define CalcLowerTriangularDet(x) calcDiagonalDet(x)

	inline long double calcCotDiagonalDet(Determinant *cdiag);
#define CalcCotUpperTriangularDet(x) calcCotDiagonalDet(x)
#define CalcCotLowerTriangularDet(x) calcCotDiagonalDet(x)

	std::vector<std::string> getPermutation(std::string *numCombination);
	unsigned int  getInverseNumber(std::string *numCombination);

	std::ostream &operator<<(std::ostream &stream, Determinant &det)
	{
		for (inSDef::Counting icount = 1; icount <= det.order(); icount++)
		{
			stream << "|";
			for (inSDef::Counting jcount = 1; jcount <= det.order(); jcount++)
			{
				stream << std::setw(4) << det(icount, jcount);
			}
			stream << " |" << std::endl;
		}
		return stream;
	}

	std::istream &operator>>(std::istream &stream, Determinant &det)
	{
		for (inSDef::Counting icount = 1; icount <= det.order(); icount++)
		{
			for (inSDef::Counting jcount = 1; jcount <= det.order(); jcount++)
			{
				std::cout << "input a" << std::setw(2) << icount << "," << std::setw(2) << jcount << "->";
				stream >> det(icount, jcount);
			}
		}
		return stream;
	}

	//long double operator*(long double k, Det &det)
	//{
	//	return k * det.calcDet();
	//}

	//long double operator*(Det &det, long double k)
	//{
	//	return k * det.calcDet();
	//}

	inline inSDef::RealNumValue calcDiagonalDet(Determinant *diag)
	{
		inSDef::Line ord = diag->order();
		inSDef::RealNumValue result = 0;
		for (inSDef::Counting orderCount = 1; orderCount <= ord; orderCount++)
		{
			result *= (*diag)(orderCount, orderCount);
		}
		return result;
	}

	inline inSDef::RealNumValue calcCotDiagonalDet(Determinant *cdiag)
	{
		inSDef::Line ord = cdiag->order();
		inSDef::RealNumValue result = 0;
		for (inSDef::Counting icount = 1, jcount = ord; icount <= ord; icount++, jcount--)
		{
			result *= (*cdiag)(icount, jcount);
		}
		result *= Determinant::negativeOnePower(ord * (ord - 1) / 2);
		return result;
	}
}

#endif