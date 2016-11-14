#include "../LAinc/DeterminantDeclaration.h"


//using does not influence client's demo in c_plus_plus file

using namespace std;
using namespace LAP;
using namespace LAP::inSDef;

Determinant::Determinant() : m_UIOrder(0), m_BChangeSymbol(false)
{

}

Determinant::Determinant(Line order) : m_UIOrder(order), m_BChangeSymbol(false)
{
	for (Counting icount = 1; icount <= order; icount++)
	{
		for (Counting jcount = 1; jcount <= order; jcount++)
		{
			this->tmpv.push_back(0);
		}
		this->m_SVecDimTwoTable.push_back(this->tmpv);
		this->tmpv.clear();
	}
	this->iter = this->m_SVecDimTwoTable.begin();
	this->it = this->iter->begin();
}

Determinant::~Determinant()
{
	
}

inline Determinant Determinant::getMinorMij(Line i, Line j)
{
	Determinant Mij(this->m_UIOrder - 1);
	Mij.iter = Mij.m_SVecDimTwoTable.begin();
	Mij.tmpv.clear();
	for (Counting row = 1; row <= this->order(); row++)
	{
		for (Counting column = 1; column <= this->order(); column++)
		{
			if (row == i || column == j)
			{
				continue;
			}
			Mij.tmpv.push_back((*this)(row, column));
		}
		if (Mij.tmpv.empty())
		{
			continue;
		}
		Mij.iter->assign(Mij.tmpv.begin(), Mij.tmpv.end());
		Mij.iter++;
		Mij.tmpv.clear();
	}
	Mij.iter = Mij.m_SVecDimTwoTable.begin();
	return Mij;
}

long double Determinant::calcMinorMij(Line i, Line j)
{
	return this->getMinorMij(i, j).calcDet();
}

long double Determinant::calcDet()
{
	RealNumValue sum = 0;
	switch (this->m_UIOrder)
	{
	case 0:
		cout << "\nA determinant can not be zero order" << endl;
		exit(EXIT_FAILURE);
		break;
	case 1:
		return (*this)(1, 1);
	case 2:
		return (*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1);
	case 3:
		return (*this)(1, 1) * (*this)(2, 2)*(*this)(3, 3) + (*this)(1, 2) * (*this)(2, 3) * (*this)(3, 1)
			+ (*this)(1, 3) * (*this)(2, 1) * (*this)(3, 2) - (*this)(1, 3) * (*this)(2, 2) * (*this)(3, 1)
			- (*this)(1, 2) * (*this)(2, 1) * (*this)(3, 3) - (*this)(1, 1) * (*this)(2, 3) * (*this)(3, 2);
	default:
		for (unsigned int i = 0; i < this->m_UIOrder; i++)
		{
			sum += negativeOnePower(1 + i) * (*this)(1, i) * this->getMinorMij(1, i).calcDet();
		}
		return sum;
	}
}

long double Determinant::calcCofactor(Line i, Line j)
{
	return negativeOnePower(i + j) * this->getMinorMij(i, j).calcDet();
}

long double Determinant::powerItself(Pow power)
{
	return pow(this->calcDet(), power);
}

void Determinant::swapIJ(Line i, Line j, bool isRow /* = true */)
{
	if (isRow)
	{
		this->iter = this->m_SVecDimTwoTable.begin();
		this->iter += i - 1;
		this->iter->swap(*(this->m_SVecDimTwoTable.begin() + j - 1));
	}
	else
	{
		//put j column into tmpv
		this->tmpv.clear();
		for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
		{
			this->tmpv.push_back((*this)(icount, j));
		}
		//replace j column with i column
		for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
		{
			(*this)(icount, j) = (*this)(icount, i);
		}
		//fill i column with tmpv
		for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
		{
			(*this)(icount, i) = this->tmpv.at(icount - 1);
		}
		this->tmpv.clear();
	}
	this->m_BChangeSymbol = !this->m_BChangeSymbol;
}

Determinant Determinant::getTranspose()
{
	Determinant dtmp(this->m_UIOrder);
	for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_UIOrder; jcount++)
		{
			dtmp(icount, jcount) = (*this)(jcount, icount);
		}
	}
	return dtmp;
}

void Determinant::transposeMyself()
{
	*this = getTranspose();
}

void Determinant::NumMultiply(const RealNumValue &real, Line line, bool isRow /* = true */)
{
	if (isRow)
	{
		for (Counting jcount = 1; jcount <= this->m_UIOrder; jcount++)
		{
			(*this)(line, jcount) *= real;
		}
	}
	else
	{
		for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
		{
			(*this)(icount, line) *= real;
		}
	}
}

void Determinant::Initialize(long double *addressOfArrayFirstElement)
{
	for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_UIOrder; jcount++)
		{
			(*this)(icount, jcount) = *(addressOfArrayFirstElement + (icount - 1) * this->m_UIOrder + jcount - 1);
		}
	}
}

void Determinant::Initialize()
{
	for (Counting icount = 1; icount <= m_UIOrder; icount++)
	{
		for (Counting jcount = 1; jcount <= m_UIOrder; jcount++)
		{
			cout << "input a" << setw(2) << icount << "," << setw(2) << jcount << "->";
			cin >> (*this)(icount, jcount);
		}
	}
}

void Determinant::displayDet()
{
	for (Counting icount = 1; icount <= m_UIOrder; icount++)
	{
		cout << "|";
		for (Counting jcount = 1; jcount <= m_UIOrder; jcount++)
		{
			cout << setw(4) << (*this)(icount, jcount);
		}
		cout << " |" << endl;
	}
}

long double Determinant::operator+()
{
	return this->calcDet();
}

long double Determinant::operator-()
{
	return -1 * this->calcDet();
}

bool Determinant::operator>(Determinant &det)
{
	return (this->calcDet() > det.calcDet());
}

bool Determinant::operator>=(Determinant &det)
{
	return (this->calcDet() >= det.calcDet());
}

bool Determinant::operator<(Determinant &det)
{
	return (this->calcDet() < det.calcDet());
}

bool Determinant::operator<=(Determinant &det)
{
	return (this->calcDet() <= det.calcDet());
}

bool Determinant::operator!=(Determinant &det)
{
	return (this->calcDet() != det.calcDet());
}

bool Determinant::operator==(Determinant &det)
{
	return (this->calcDet() == det.calcDet());
}

long double Determinant::operator+(Determinant &det)
{
	return this->calcDet() + det.calcDet();
}

long double Determinant::operator-(Determinant &det)
{
	return this->calcDet() - det.calcDet();
}

long double Determinant::operator*(Determinant &det)
{
	return this->calcDet() * det.calcDet();
}

long double Determinant::operator/(Determinant &det)
{
	RealNumValue dive = 0;
	dive = det.calcDet();
	if (dive)
	{
		cout << "\nError! divisor is zero, so result is 0 to indicate failure..." << endl;
		return 0;
	}
	return this->calcDet() / dive;
}

Determinant & Determinant::operator=(Determinant &det)
{
	if (this->order() != det.order())
	{
		cout << "\nError! assign elements with a determinant whose order is not same." << endl;
		cout << "so return the left determinant to indicate failure..." << endl;
		return *this;
	}
	for (Counting icount = 1; icount <= this->m_UIOrder; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_UIOrder; jcount++)
		{
			(*this)(icount, jcount) = det(icount, jcount);
		}
	}
	this->m_BChangeSymbol = det.m_BChangeSymbol;
	this->tmpv.clear();
	this->iter = this->m_SVecDimTwoTable.begin();
	return *this;
}