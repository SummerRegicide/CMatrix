#ifndef REAL_NUMBER_SYS_MATRIX_FRAME_H
#define REAL_NUMBER_SYS_MATRIX_FRAME_H

#include "RealSMatrixDeclaration.h"

namespace LAP
{
	RealSysMatrix createIdentityMatrix(inSDef::Line order);
	RealSysMatrix createCotIdentityMatrix(inSDef::Line order);
	RealSysMatrix createDiagonalMatrix(inSDef::Line order, inSDef::RealNumValue lambda);
	RealSysMatrix createCotDiagonalMatrix(inSDef::Line order, inSDef::RealNumValue lambda);
	
	std::ostream &operator<<(std::ostream &stream, RealSysMatrix &matrix);
	std::istream &operator>>(std::istream &stream, RealSysMatrix &matrix);
	
	RealSysMatrix operator*(RealSysMatrix &matrix, long double &k);
	RealSysMatrix operator*(long double &k, RealSysMatrix &matrix);
	RealSysMatrix operator+(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
	RealSysMatrix operator-(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
	RealSysMatrix operator*(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
	RealSysMatrix &operator+=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
	RealSysMatrix &operator-=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
	RealSysMatrix &operator*=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);

}

#endif