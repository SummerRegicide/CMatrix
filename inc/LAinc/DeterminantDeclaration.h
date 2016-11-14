#ifndef _DETDECLARATION_H__
#define _DETDECLARATION_H__

#include "GeneralDef.h"

namespace LAP
{
	class Determinant
	{
	public:

		//Though it's a no-argument constructor, do not use
		Determinant();
		Determinant(inSDef::Line order);
		//do not use copy constructor
		/*Determinant(Determinant &);*/

		~Determinant();

		//An n*n array is okay
		void Initialize(long double *addressOfArrayFirstElement);

		//input by yourself
		void Initialize();

		//print formatted
		void displayDet();

		unsigned int order() const
		{
			return m_UIOrder;
		}

		unsigned int sizeOfElements() const
		{
			return m_UIOrder * m_UIOrder;
		}

		inSDef::RealNumValue calcDet();

		Determinant getTranspose();
		inline Determinant getMinorMij(inSDef::Line i, inSDef::Line j);
		inSDef::RealNumValue calcMinorMij(inSDef::Line i, inSDef::Line j);
		inSDef::RealNumValue calcCofactor(inSDef::Line i, inSDef::Line j);
		void transposeMyself();
		void swapIJ(inSDef::Line i, inSDef::Line j, bool isRow = true);

		void NumMultiply(const inSDef::RealNumValue &real, inSDef::Line line, bool isRow = true);
		inSDef::RealNumValue powerItself(inSDef::Pow power);

		//Overload operator functions for output and input stream
		inline friend std::ostream &operator<<(std::ostream &stream, Determinant &det);
		inline friend std::istream &operator>>(std::istream &stream, Determinant &det);

		//Overload operator functions for unary + and -
		inSDef::RealNumValue operator+();
		inSDef::RealNumValue operator-();

		//Overload operator functions for relational operators
		bool operator>(Determinant &det);
		bool operator<(Determinant &det);
		bool operator>=(Determinant &det);
		bool operator==(Determinant &det);
		bool operator<=(Determinant &det);
		bool operator!=(Determinant &det);

		//Overload operator function for conversion
		operator long double()
		{
			return this->calcDet();
		}

		//Overload operator function for assignment
		//It's unsafe. Avoid using.
		Determinant &operator=(Determinant &det);

		//Overload operator function for multiplying by a number
		//They are repetitive.
		//friend inSDef::RealNumValue operator*(long double k, Determinant &det);
		//friend inSDef::RealNumValue operator*(Determinant &det, long double k);

		//Overload operator functions for arithmetic operators among determinants
		//It's safe to divide.Do not worry.
		inSDef::RealNumValue operator+(Determinant &det);
		inSDef::RealNumValue operator-(Determinant &det);
		inSDef::RealNumValue operator*(Determinant &det);
		inSDef::RealNumValue operator/(Determinant &det);

		//Overload operator functions for accessing A(i,j)
		inSDef::RealNumValue & operator()(inSDef::Line i, inSDef::Line j)
		{
			if (i == 0 || j == 0)
			{
				exit(EXIT_FAILURE);
			}
			if (i > m_UIOrder || j > m_UIOrder)
			{
				exit(EXIT_FAILURE);
			}
			this->iter = this->m_SVecDimTwoTable.begin();
			this->iter += i - 1;
			this->it = this->iter->begin();
			this->it += j - 1;
			return *it;
		}

		inline friend inSDef::RealNumValue calcCotDiagonalDet(Determinant *cdiag);

	/*protected:*/
	private:
		inSDef::Symbol m_BChangeSymbol;
		inSDef::Line m_UIOrder;
		inSDef::VecSetIter iter;
		inSDef::VecIter it;
		inSDef::Vec tmpv;
		inSDef::VecSet m_SVecDimTwoTable;
		static inSDef::PosOrNeg negativeOnePower(inSDef::Pow power)
		{
			return (power % 2) ? -1 : 1;
		}
		inSDef::RealNumValue getElementIJ(inSDef::Line i, inSDef::Line j)
		{
			return (*this)(i, j);
		}
			
	};

	typedef Determinant Det;
}

#endif