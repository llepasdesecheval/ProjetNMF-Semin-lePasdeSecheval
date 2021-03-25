//
//  Matrix.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 25/03/2021.
//

#include "Matrix.hpp"

#include <functional>

Matrix::Matrix(std::size_t nRows, std::size_t nCols, bool zeros) : m_nRows(nRows), m_nCols(nCols)
{
    if (zeros)
        m_data = std::vector<double>(m_nRows * m_nCols, 0);
    else
        m_data = std::vector<double>(m_nRows * m_nCols);
}

Matrix::Matrix(std::size_t nRows, std::size_t nCols, double value) : m_nRows(nRows), m_nCols(nCols), m_data(m_nRows * m_nCols, value) {}

Matrix::Matrix(std::string data)
{
    char separator;
    // read once to compute the size
    m_nRows = 1;
    m_nCols = 1;
    std::size_t n = data.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (data[i] == ';')
        {
            // We add a new row when encountering ';'
            m_nRows++;
        }
        if (m_nRows==1 && data[i] == ',')
        {
            // As long as we did not meet ';', the number of columns increases
            m_nCols++;
        }
    }
    
    // now check we can read the string
    std::stringstream ss1;
    ss1.str(data);
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        for (std::size_t j = 0; j < m_nCols; ++j)
        {
            double ignored;
            ss1 >> ignored;
            ss1 >> separator;
            if ((j == nCols() - 1) && (i < nRows() - 1))
            {
                ASSERT(separator == ';');
            }
            else if (j < nCols() - 1)
            {
                ASSERT(separator == ',');
            }
        }
    }
    
    // allocate memory now we know nothing will go wrong
    std::stringstream ss;
    ss.str(data);
    m_data = std::vector<double>(m_nRows * m_nCols);
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        for (std::size_t j = 0; j < m_nCols; ++j)
        {
            ss >> m_data[offset(i, j)];
            ss >> separator;
        }
    }
}

Matrix::Matrix(std::size_t nRows, std::size_t nCols, std::vector<double> data)
{
    ASSERT(data.size() == nRows * nCols);
    m_nRows = nRows;
    m_nCols = nCols;
    m_data = data;
}

Matrix Matrix::Identity(std::size_t dim)
{
    Matrix res(dim, dim);
    for (std::size_t i = 0; i < dim; ++i)
    {
        res.m_data[i * (dim + 1)] = 1;
    }
    return res;
}

Matrix Matrix::transpose() const
{
    std::size_t m = m_nRows;
    std::size_t n = m_nCols;
    std::vector<double> temp(m * n);
    
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            temp[j * m + i] = m_data[offset(i, j)];
        }
    }
    Matrix res(n, m, temp);
    return res;
}

Matrix zeros(std::size_t nRows, std::size_t nCols)
{
    return Matrix(nRows, nCols);
}

Matrix ones(std::size_t nRows, std::size_t nCols)
{
    return Matrix(nRows, nCols, 1.0);
}

Matrix Cholesky(const Matrix& m)
{
    std::size_t dim = m.nRows();
    ASSERT(dim == m.nCols()); // Cholesky works only for square matrices
    Matrix res(dim, dim, true);
    
    for (std::size_t i = 0; i < dim; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            double s = m(i, j);
            for (std::size_t k = 0; k < j; ++k)
            {
                s -= res(i, k) * res(j, k);
            }
            res(i, j) = s / res(j, j);
        }
        double s = m(i, i);
        for (std::size_t k = 0; k < i; ++k)
        {
            s -= res(i, k) * res(i, k);
        }
        // m must be positive-definite
        ASSERT(s >= 0);
        res(i, i) = sqrt(s);
    }
    
    return res;
}

void printMatrix(const Matrix& m)
{
    std::size_t limi = m.nRows() - 1;
    std::size_t limj = m.nCols();
    std::cout << "[";
    for (std::size_t i = 0; i < limi; ++i)
    {
        std::cout << "[";
        for (std::size_t j = 0; j < limj; ++j)
        {
            printf("% 8g", m(i, j));
        }
        std::cout << "]\n ";
    }
    std::cout << "[";
    for (std::size_t j = 0; j < limj; ++j)
    {
        printf("% 8g", m(limi, j));
    }
    std::cout << "]]" << std::endl;
}

Matrix& Matrix::operator+=(const Matrix& rhs)
{
    std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(), std::plus<double>());
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
    std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(), std::minus<double>());
    return *this;
}

Matrix& Matrix::operator+=(double rhs)
{
    std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double arg) { return arg + rhs; });
    return *this;
}

Matrix& Matrix::operator-=(double rhs)
{
    std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double arg) { return arg - rhs; });
    return *this;
}

Matrix& Matrix::operator*=(double rhs)
{
    std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double arg) { return arg * rhs; });
    return *this;
}

Matrix& Matrix::operator/=(double rhs)
{
    ASSERT(rhs != 0);
    std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double arg) { return arg / rhs; });
    return *this;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
    ASSERT((lhs.nRows() == rhs.nRows()) && (lhs.nCols() == rhs.nCols()));
    Matrix res(lhs);
    res += rhs;
    
    return res;
}

Matrix operator+(const Matrix& lhs, double rhs)
{
    Matrix res(lhs);
    res += rhs;
    
    return res;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
    ASSERT(lhs.nRows() == rhs.nRows() && lhs.nCols() == rhs.nCols());
    Matrix res(lhs);
    res -= lhs;
    
    return res;
}

Matrix operator-(const Matrix& lhs, double rhs)
{
    Matrix res(lhs);
    res -= rhs;
    
    return res;
}

Matrix operator-(double lhs, const Matrix& rhs)
{
    Matrix res(rhs);
    res *= -1.0;
    res += rhs;
    
    return res;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
    std::size_t m = lhs.nRows();
    std::size_t r = lhs.nCols();
    // The number of columns of the left matrix must equal the number of lines of the right matrix
    ASSERT(rhs.nRows() == r);
    std::size_t n = rhs.nCols();
    Matrix res(m, n, true);
    
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            // k running through a's lines & b's columns
            for (std::size_t k = 0; k < r; ++k)
            {
                res(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return res;
}

Matrix operator*(const Matrix& lhs, double rhs)
{
    Matrix res(lhs);
    res *= rhs;
    return res;
}

std::vector<double> operator*(const Matrix& lhs, const std::vector<double> rhs)
{
    std::size_t m = lhs.nRows();
    std::size_t n = lhs.nCols();
    ASSERT(rhs.size() == n);
    std::vector<double> res(m, 0);
    
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            res[i] += lhs(i, j) * rhs[j];
        }
    }
    return res;
}

std::vector<double> operator*(const std::vector<double> lhs, const Matrix& rhs)
{
    std::size_t m = rhs.nRows();
    std::size_t n = rhs.nCols();
    ASSERT(lhs.size() == m);
    std::vector<double> res(n, 0);
    
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            res[j] += lhs[i] * rhs(i, j);
        }
    }
    return res;
}

Matrix Matrix::Cofactor(std::size_t i, std::size_t j) const
{
    ASSERT(m_nRows == m_nCols);
    std::size_t k = 0;
    std::size_t l = 0;
    Matrix res(m_nRows - 1, m_nCols - 1);
    
    // Looping for each element of the matrix
    for (std::size_t row = 0; row < m_nRows; ++row)
    {
        for (std::size_t col = 0; col < m_nCols; ++col)
        {
            // Copying into temporary matrix only those element which are not in given row and column
            if (row != i && col != j)
            {
                res(k, l++) = m_data[offset(row, col)];
                // Row is filled, so increase row index and reset col index
                if (l == m_nCols - 1)
                {
                    l = 0;
                    ++k;
                }
            }
        }
    }
    
    return res;
}

double Matrix::Determinant() const
{
    ASSERT(m_nRows == m_nCols);
    double res = 0;
    
    //  Base case : if matrix contains single element
    if (m_nRows == 1)
        return m_data[0];
    
    Matrix temp(m_nRows - 1, m_nCols - 1); // To store cofactors
    double sign = 1.0;  // To store sign multiplier
    
    // Iterate for each element of first row
    for (std::size_t j = 0; j < m_nCols; ++j)
    {
        // Getting Cofactor of fth element
        temp = Cofactor(0, j);
        res += sign * m_data[j] * temp.Determinant();
        // Sign switches at each step
        sign *= -1;
    }
    
    return res;
}

Matrix Matrix::Adjugate() const
{
    ASSERT(m_nRows == m_nCols);
    Matrix res(m_nRows, m_nCols);
    
    if (m_nRows == 1)
    {
        res.m_data[0] = 1;
    }
    
    int sign = 1;
    // Used to store cofactors
    Matrix temp(m_nRows - 1, m_nCols - 1);
    
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        for (std::size_t j = 0; j < m_nCols; ++j)
        {
            temp = Cofactor(i, j);
            // sign of adj[j][i] positive if sum of row and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            
            // Interchanging rows and columns to get the transpose of the cofactor matrix
            res(j, i) = sign * temp.Determinant();
        }
    }
    
    return res;
}

Matrix Matrix::Inverse() const
{
    ASSERT(m_nRows == m_nCols);
    double det = Determinant();
    ASSERT(det != 0);
    Matrix adj = Adjugate();
    Matrix res(m_nRows, m_nCols, false);
    
    // A^{-1} = Adjugate(A)/Determinant(A)
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        for (std::size_t j = 0; j < m_nCols; ++j)
        {
            res(i, j) = adj(i, j)/det;
        }
    }
    
    return res;
}

void Matrix::assertEquals(const Matrix& other, double tolerance) const
{
    ASSERT((other.nRows() == nRows()) && (other.nCols() == nCols()));
    for (std::size_t i = 0; i < nRows(); ++i)
    {
        for (std::size_t j = 0; j < nCols(); ++j)
        {
            double expected = m_data[offset(i, j)];
            double actual = other(i, j);
            if (fabs(expected - actual) > tolerance)
            {
                std::stringstream s;
                s << "ASSERTION FAILED\n";
                s << "Mismatch at cell (" << i << ", " << j <<":\n";
                s << "Expected: " << expected << "; Actual: " << actual << ".\n";
                printMatrix(*this);
                printMatrix(other);
                // Prints where the error was thrown
                INFO(s.str());
                // Constructs an exception with explanatory string as argument
                throw std::runtime_error(s.str());
            }
        }
    }
}

std::vector<double> Matrix::ExtractRow(std::size_t i) const
{
    std::vector<double> res(m_nCols);
    std::copy(m_data.begin() + i * m_nCols, m_data.begin() + (i + 1) * m_nCols, res.begin());
    return res;
}

void Matrix::setRow(std::size_t i, const std::vector<double>& data)
{
    ASSERT(data.size() == m_nCols);
    std::copy(data.begin(), data.end(), m_data.begin() + i * m_nCols);
    
}

std::vector<double> Matrix::asRow() const
{
    ASSERT(m_nRows == 1);
    return m_data;
}

std::vector<double> Matrix::ExtractColumn(std::size_t j) const
{
    std::vector<double> res(m_nRows);
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        res[i] = m_data[offset(i, j)];
    }
    return res;
}

void Matrix::setColumn(std::size_t j, const std::vector<double>& data)
{
    ASSERT(data.size() == m_nRows);
    for (std::size_t i = 0; i < m_nRows; ++i)
    {
        m_data[offset(i, j)] = data[i];
    }
}

std::vector<double> Matrix::asColumn() const
{
    ASSERT(m_nCols == 1);
    return m_data;
}

static void testConstructor()
{
    Matrix m(3,8);
    ASSERT(m.nRows() == 3);
    ASSERT(m.nCols() == 8);
    for (std::size_t i = 0; i < m.nRows(); ++i)
    {
        for (std::size_t j = 0; j < m.nCols(); ++j)
        {
            ASSERT(m(i, j) == 0.0);

            m(i,j) = i + j;
            ASSERT(m(i,j) == i+j);
        }
    }
    
    Matrix m1(4, 4, 1.0);
    Matrix m2 = ones(4, 4);
    m1.assertEquals(m2, 1E-6);
}

static void testIdentity()
{
    Matrix I = Matrix::Identity(3);
    printMatrix(I);
}

static void testTranspose()
{
    Matrix m("3;2;1");
    Matrix expected("3,2,1");
    expected.assertEquals(m.transpose(), 1E-6);
}

static void testDeterminant()
{
    Matrix m1("2,3;1,2");
    ASSERT(m1.Determinant() == 1);
    
    Matrix m2("1,3,5;3,4,8;5,9,13");
    ASSERT(m2.Determinant() == 18);
    
    Matrix m3("4,15,20,15;20,10,13,9;19,7,6,2;14,10,2,14");
    ASSERT(m3.Determinant() == 15228);
}

static void testAdditionAndSubtrationOperators()
{
    Matrix z = zeros(3, 2);
    Matrix m = ones(3, 2);
    Matrix n = ones(3, 2);

    m.assertEquals(n, 1E-6);
    (1 + m).assertEquals(2 * n, 1E-6);
    (1 + m).assertEquals(n * 2, 1E-6);
    (m + 1).assertEquals(n * 2, 1E-6);
    Matrix temp = m + m + n;
    //temp += n;
    temp.assertEquals(n * 3, 1E-6);
    (m - 1).assertEquals(z, 1E-6);
    (1 - m).assertEquals(z, 1E-6);
    (m - m).assertEquals(z, 1E-6);
}

static void testAssignmentOperators()
{
    const Matrix u = ones(3,2);
    Matrix m = ones(3,2);
    
    m = ones(3,2);
    m += 1;
    m.assertEquals(u + u, 0.001);

    m += 2 * u;
    m.assertEquals(4 * u, 0.001);

    m -= 2 * u;
    m.assertEquals(2 * u, 0.001);

    m -= 1;
    m.assertEquals(u, 0.001);

    m *= 8;
    m.assertEquals(8 * u, 0.001);
}

static void testReadFromString()
{
    Matrix m("1,2,3;4,5,6");
    ASSERT(m.nRows() == 2);
    ASSERT(m.nCols() == 3);
    double count = 1.0;
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            ASSERT_APPROX_EQUAL(m(i, j), count++, 1E-6);
        }
    }
}

static void testMatrixMultiplication()
{
    Matrix a("1,2,3;4,5,6");
    Matrix b("1,2;3,4;5,6");
    Matrix product = a*b;
    Matrix expected("22,28;49,64");
    product.assertEquals(expected, 1E-6);
}

static void testMatrixVectorMultiplication()
{
    std::vector<double> v(4);
    Matrix I = Matrix::Identity(4);
    std::vector<double> v1 = v * I;
    std::vector<double> v2 = I * v;
    ASSERT(v1 == v);
    ASSERT(v2 == v);
    
    Matrix zm(4, 4);
    std::vector<double> zv(4, 0);
    std::vector<double> v3 = zm * v;
    std::vector<double> v4 = v * zm;
    ASSERT(v3 == zv);
    ASSERT(v4 == zv);
    
    Matrix m("1,2,3;4,5,6");
    std::vector<double> ov1(3, 1);
    std::vector<double> ov2(2, 1);
    std::vector<double> r1 = {6, 15};
    std::vector<double> r2 = {5, 7, 9};
    std::vector<double> v5 = m * ov1;
    std::vector<double> v6 = ov2 * m;
    ASSERT(v5 == r1);
    ASSERT(v6 == r2);
}

static void testCholesky()
{
    Matrix m("3,1,2;1,4,-1;2,-1,5");
    Matrix c = Cholesky(m);
    Matrix product = c * c.transpose();
    m.assertEquals(product, 1E-6);
}

static void testInverse()
{
    Matrix m("3,1,2;1,4,-1;2,-1,5");
    Matrix product1 = m * m.Inverse();
    Matrix product2 = m.Inverse() * m;
    Matrix I = Matrix::Identity(m.nRows());
    product1.assertEquals(I, 1E-6);
    product2.assertEquals(I, 1E-6);
}

static void testExtractors()
{
    Matrix m("1,2,3;4,5,6;7,8,9");
    std::vector<double> row2 = m.ExtractRow(1);
    std::vector<double> expectedRow = {4, 5, 6};
    
    std::vector<double> col3 = m.ExtractColumn(2);
    std::vector<double> expectedCol = {3, 6, 9};
    
    for (std::size_t i = 0; i < 3; ++i)
    {
        ASSERT(row2[i] == expectedRow[i]);
        ASSERT(col3[i] == expectedCol[i]);
    }
}

static void testRowAndColSetters()
{
    Matrix m("1,2,3,4;5,6,7,8;9,10,11,12");
    std::vector<double> row2 = {13, 14, 15, 16};
    m.setRow(1, row2);
    std::vector<double> newRow2 = m.ExtractRow(1);
    ASSERT(row2 == newRow2);
    
    std::vector<double> col3 = {10, 15, 20};
    m.setColumn(2, col3);
    std::vector<double> newCol3 = m.ExtractColumn(2);
    ASSERT(col3 == newCol3);
}

void testMatrix()
{
    TEST(testConstructor);
    TEST(testIdentity);
    TEST(testTranspose);
    TEST(testDeterminant);
    TEST(testInverse);
    TEST(testAdditionAndSubtrationOperators);
    TEST(testAssignmentOperators);
    TEST(testReadFromString);
    TEST(testMatrixMultiplication);
    TEST(testMatrixVectorMultiplication);
    TEST(testCholesky);
    TEST(testExtractors);
    TEST(testRowAndColSetters);
}
