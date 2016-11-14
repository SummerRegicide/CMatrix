#include "../LAinc/RealSMatrixDeclaration.h"
#include <algorithm>

/*
some flags:
	//defective		have not finished
	//slow			not efficient
	//unsafe		need improving
*/
using namespace std;
using namespace LAP;
using namespace LAP::inSDef;

template <typename _Type>
class MultiVal
{
public:
	MultiVal(const _Type &val) : factor(val) {}
	void operator()(_Type &ele)
	{
		ele *= factor;
	}
private:
	_Type factor;
};

template <typename _Type>
class PlusVal
{
public:
	PlusVal(const _Type &val) : factor(val) {}
	void operator()(_Type &ele)
	{
		ele += factor;
	}
private:
	_Type factor;
};

RealSysMatrix::RealSysMatrix() : m_CouRow(0), m_CouColumn(0)
{
	//defective
	cout << "\nError! no-argument constructor is not allowed..." << endl;
	exit(EXIT_FAILURE);
}

RealSysMatrix::RealSysMatrix(Line m, Line n) : m_CouRow(m), m_CouColumn(n)
{
	//defective
	if (!m)
	{
		cout << "\nError! 0 row is not allowed..." << endl;
		exit(EXIT_FAILURE);
	}
	if (!n)
	{
		cout << "\nError! 0 column is not allowed..." << endl;
		exit(EXIT_FAILURE);
	}
	this->tmpv.clear();
	this->m_SVecDimTwoTable.clear();
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= m_CouColumn; jcount++)
		{
			this->tmpv.push_back(0);
		}
		this->m_SVecDimTwoTable.push_back(this->tmpv);
		this->tmpv.clear();
	}
	this->m_BIsSquare = (m == n);
	this->m_BIsSymmetric = false;
	this->m_BIsInvertible = false;
}

RealSysMatrix::RealSysMatrix(const RealSysMatrix &matrix) : m_CouRow(matrix.Row()), m_CouColumn(matrix.Column())
{
	this->m_BIsSquare = matrix.m_BIsSquare;
	this->m_BIsSymmetric = matrix.m_BIsSymmetric;
	this->m_BIsInvertible = matrix.m_BIsInvertible;
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) = matrix(icount, jcount);
		}
	}
	this->tmpv.clear();
	this->UpToDate();
}

unsigned int RealSysMatrix::SizeOfElement() const
{
	return this->m_CouRow * this->m_CouColumn;
}

RealSysMatrix::~RealSysMatrix()
{
	this->tmpv.clear();
	this->m_SVecDimTwoTable.clear();
}

Line RealSysMatrix::Row() const
{
	return this->m_CouRow;
}

Line RealSysMatrix::Column() const
{
	return this->m_CouColumn;
}

void RealSysMatrix::Initialize(long double *addressOfArrayFirstElement)
{
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) = *(addressOfArrayFirstElement + (icount - 1) * this->m_CouRow + jcount - 1);
		}
	}
	this->UpToDate();
}

void RealSysMatrix::Initialize(long double **arr)
{
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) = *(*arr);
			++arr;
		}
	}
	this->UpToDate();
}

void RealSysMatrix::Initialize()
{
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			cout << "input a" << setw(2) << icount << "," << setw(2) << jcount << "->";
			cin >> (*this)(icount, jcount);
		}
	}
	this->UpToDate();
}

bool RealSysMatrix::isSquare()
{
	return (this->m_CouRow == this->m_CouColumn);
}

bool RealSysMatrix::isSymmetric()
{
	//defective
	if (!this->m_BIsSquare)
	{
		return false;
	}
	return ((*this) == this->getTranspose());
}

bool RealSysMatrix::isInvertible()
{
	//defective
	if (!this->m_BIsSquare)
	{
		return false;
	}
	return (this->Rank() == this->getLess());
}

unsigned int RealSysMatrix::Rank()
{
	return 0;
}

long double RealSysMatrix::Trace()
{
	//defective
	if (!this->isSquare())
	{
		cout << "\nError! matrix is not square." << endl;
		cout << "return 0 to indicate failure..." << endl;
		return 0;
	}
	long double sum = 0;
	for (Counting ct = 1; ct <= this->m_CouRow; ct++)
	{
		sum += (*this)(ct, ct);
	}
	return sum;
}

RealNumValue & RealSysMatrix::operator()(Line i, Line j)
{
	//defective
	if (!i)
	{
		cout << "\nError! can not access 0 row..." << endl;
		exit(EXIT_FAILURE);
	}
	if (!j)
	{
		cout << "\nError! can not access 0 column..." << endl;
		exit(EXIT_FAILURE);
	}
	/*VecSetIter iter = this->m_SVecDimTwoTable.begin();
	iter += i - 1;
	VecIter it = iter->begin();
	it += j - 1;
	return *it;*/
	return this->m_SVecDimTwoTable.at(i - 1).at(j - 1);
}

const RealNumValue & RealSysMatrix::operator()(Line i, Line j) const
{
	//defective
	if (!i)
	{
		cout << "\nError! can not access 0 row..." << endl;
		exit(EXIT_FAILURE);
	}
	if (!j)
	{
		cout << "\nError! can not access 0 column..." << endl;
		exit(EXIT_FAILURE);
	}
	cVecSetIter iter = this->m_SVecDimTwoTable.begin();
	iter += i - 1;
	cVecIter it = iter->begin();
	it += j - 1;
	return *it;
}

Det RealSysMatrix::getDet()
{
	//defective
	if (!this->m_BIsSquare)
	{

	}
	Determinant tmp_det(this->m_CouRow);
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			tmp_det(icount, jcount) = (*this)(icount, jcount);
		}
	}
	return tmp_det;
}

RealSysMatrix RealSysMatrix::getTranspose()
{
	RealSysMatrix tmp_trans(this->m_CouColumn, this->m_CouRow);
	for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
	{
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			tmp_trans(jcount, icount) = (*this)(icount, jcount);
		}
	}
	return tmp_trans;
}

void RealSysMatrix::Transpose()
{
	RealSysMatrix tmp_matrix(*this);
	this->m_SVecDimTwoTable.clear();
	for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
	{
		this->tmpv.clear();
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			this->tmpv.push_back(tmp_matrix(jcount, icount));
		}
		this->m_SVecDimTwoTable.push_back(this->tmpv);
	}
	this->tmpv.clear();
	this->UpToDate();
}

RealSysMatrix RealSysMatrix::getAdjugate()
{
	//defective
	if (!this->m_BIsSquare)
	{
	}
	RealSysMatrix tmp_matrix(this->m_CouRow, this->m_CouColumn);
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			tmp_matrix(jcount, icount) = this->getDet().calcCofactor(icount, jcount);
		}
	}
	return tmp_matrix;
}

void RealSysMatrix::Adjugate()
{
	*this = this->getAdjugate();
}

RealSysMatrix RealSysMatrix::getInverse()
{
	RealNumValue det = this->getDet().calcDet();
	RealSysMatrix tmp_matrix = this->getAdjugate();
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			tmp_matrix(icount, jcount) *= det;
		}
	}
	return tmp_matrix;
}

void RealSysMatrix::Inverse()
{
	RealNumValue det = this->getDet().calcDet();
	this->Adjugate();
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) *= det;
		}
	}
}

void RealSysMatrix::SwapIJ(inSDef::Line i, inSDef::Line j, bool isRow /* = true */)
{
	if (isRow)
	{
		VecSetIter outIter = this->m_SVecDimTwoTable.begin();
		outIter += i - 1;
		outIter->swap(*(this->m_SVecDimTwoTable.begin() + j - 1));
	} 
	else
	{
		this->tmpv.clear();
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			this->tmpv.push_back((*this)(icount, j));
		}
		//replace j column with i column
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			(*this)(icount, j) = (*this)(icount, i);
		}
		//fill i column with tmpv
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			(*this)(icount, i) = this->tmpv.at(icount - 1);
		}
		this->tmpv.clear();
	}
}

void RealSysMatrix::KMultiply(const long double &k)
{
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) *= k;
		}
	}
}

void RealSysMatrix::KTimesPlus(const long double &k, inSDef::Line i, inSDef::Line j, bool isRow /* = true */)
{
	//j = k * i + j
	if (isRow)
	{
		VecSetIter outIter = this->m_SVecDimTwoTable.begin() + i - 1;
		this->tmpv.assign(outIter->begin(), outIter->end());
		for_each(this->tmpv.begin(), this->tmpv.end(), MultiVal<RealNumValue>(k));
		outIter = this->m_SVecDimTwoTable.begin() + j - 1;
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			outIter->at(jcount - 1) += this->tmpv.at(jcount - 1);
		}
		this->tmpv.clear();
	} 
	else
	{
		this->tmpv.clear();
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			this->tmpv.push_back(this->m_SVecDimTwoTable.at(icount - 1).at(i - 1) * k);
		}
		for (Counting icount = 1; icount <= this->m_CouRow; icount++)
		{
			this->m_SVecDimTwoTable.at(icount - 1).at(j - 1) += this->tmpv.at(icount - 1);
		}
		this->tmpv.clear();
	}
}

RealSysMatrix RealSysMatrix::operator+()
{
	return *this;
}

RealSysMatrix RealSysMatrix::operator-()
{
	RealSysMatrix tmp_matrix(*this);
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			tmp_matrix(icount, jcount) *= -1;
		}
	}
	return tmp_matrix;
}

bool RealSysMatrix::operator==(RealSysMatrix &matrix)
{
	//defective
	if (this->m_CouRow != matrix.m_CouRow)
	{
		cout << "\nError! they have different rows..." << endl;
		exit(EXIT_FAILURE);
	}
	if (this->m_CouColumn != matrix.m_CouColumn)
	{
		cout << "\nError! they have different columns..." << endl;
		exit(EXIT_FAILURE);
	}
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			if ((*this)(icount, jcount) != matrix(icount, jcount))
			{
				return false;
			}
		}
	}
	return true;
}

bool RealSysMatrix::operator!=(RealSysMatrix &matrix)
{
	return !(*this == matrix);
}

RealSysMatrix RealSysMatrix::operator!()
{
	return this->getInverse();
}

RealSysMatrix &RealSysMatrix::operator=(const RealSysMatrix &matrix)
{
	//defective
	if (this->m_CouRow != matrix.m_CouRow)
	{
		cout << "\nError! they have different rows..." << endl;
		exit(EXIT_FAILURE);
	}
	if (this->m_CouColumn != matrix.m_CouColumn)
	{
		cout << "\nError! they have different columns..." << endl;
		exit(EXIT_FAILURE);
	}
	for (Counting icount = 1; icount <= this->m_CouRow; icount++)
	{
		for (Counting jcount = 1; jcount <= this->m_CouColumn; jcount++)
		{
			(*this)(icount, jcount) = matrix(icount, jcount);
		}
	}
	this->UpToDate();
	return *this;
}