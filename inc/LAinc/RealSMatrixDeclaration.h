#ifndef REALNUMBERSYSTEM_MATRIX_DEC_H__
#define REALNUMBERSYSTEM_MATRIX_DEC_H__

#include "GeneralDef.h"
#include "Determinant.h"

//#pragma comment(lib, "../LAlib/lapdet.lib")

namespace LAP
{
	class RealSysMatrix
	{
	public:
		RealSysMatrix();
		RealSysMatrix(inSDef::Line m, inSDef::Line n);
		RealSysMatrix(const RealSysMatrix &matrix);
		~RealSysMatrix();

		void Initialize();
		void Initialize(long double *addressOfArrayFirstElement);
		void Initialize(long double **arr);
		
		inSDef::Line Row() const;
		inSDef::Line Column() const;
		unsigned int SizeOfElement() const;

		unsigned int Rank();
		bool isSquare();
		long double Trace();
		bool isInvertible();
		bool isSymmetric();
		void displayMatrix();

		Det getDet();
		RealSysMatrix getTranspose();
		void Transpose();
		RealSysMatrix getAdjugate();
		void Adjugate();
		RealSysMatrix getInverse();
		void Inverse();
		RealSysMatrix getPower(inSDef::Pow power);
		void Power(inSDef::Pow power);

		void SwapIJ(inSDef::Line i, inSDef::Line j, bool isRow = true);
		void KMultiply(const long double &k);
		void KTimesPlus(const long double &k, inSDef::Line i, inSDef::Line j, bool isRow = true);

		RealSysMatrix operator+();
		RealSysMatrix operator-();

		friend std::ostream &operator<<(std::ostream &stream, RealSysMatrix &matrix);
		friend std::istream &operator>>(std::istream &stream, RealSysMatrix &matrix);

		bool operator==(RealSysMatrix &matrix);
		bool operator!=(RealSysMatrix &matrix);

		RealSysMatrix &operator=(const RealSysMatrix &matrix);

		friend RealSysMatrix operator*(RealSysMatrix &matrix, long double &k);
		friend RealSysMatrix operator*(long double &k, RealSysMatrix &matrix);
		friend RealSysMatrix operator+(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
		friend RealSysMatrix operator-(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
		friend RealSysMatrix operator*(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
		friend RealSysMatrix &operator+=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
		friend RealSysMatrix &operator-=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);
		friend RealSysMatrix &operator*=(RealSysMatrix &lmatrix, RealSysMatrix &rmatrix);

		RealSysMatrix operator!();

		inSDef::RealNumValue &operator()(inSDef::Line i, inSDef::Line j);
		const inSDef::RealNumValue &operator()(inSDef::Line i, inSDef::Line j) const;

	protected:
		bool m_BIsSquare;
		bool m_BIsInvertible;
		inSDef::Counting m_UIRank;
		const inSDef::Line m_CouRow;		//do not try to change the size
		const inSDef::Line m_CouColumn;		//so row and column are const.
		inSDef::RealNumValue m_RNSysTrace;
		//inSDef::VecSetIter outIter;
		//inSDef::VecIter inIter;
		inSDef::Vec tmpv;
		inSDef::VecSet m_SVecDimTwoTable;

	private:
		bool m_BIsSymmetric;
		void UpToDate()
		{
			this->m_BIsSquare = this->isSquare();
			this->m_BIsSymmetric = this->isSymmetric();
			this->m_BIsInvertible = this->isInvertible();
		}
		inSDef::Line getLess()
		{
			return (this->m_CouRow >= this->m_CouColumn) ? m_CouColumn : m_CouRow;
		}
	};

}

#endif