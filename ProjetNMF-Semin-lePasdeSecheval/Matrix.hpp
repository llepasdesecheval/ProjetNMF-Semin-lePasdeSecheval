//
//  Matrix.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 25/03/2021.
//

#if !defined MATRIX_HPP
#define MATRIX_HPP

#include "stdafx.hpp"

class Matrix
{
public:
    /**
     Constructor; all elements initialised at 0
     @param nRows Number of rows
     @param nCols Number of columns
     @param zeros If true, initialises the matrix at 0; default is true
     */
    Matrix(std::size_t nRows, std::size_t nCols, bool zeros = true);
    
    /**
     Constructor initialising all elements to some value
     @param nRows Number of rows
     @param nCols Number of columns
     @param value Initialisation value
     */
    Matrix(std::size_t nRows, std::size_t nCols, double value);
    
    /**
     Constructor from a std::string
     @param data String containing the matrix
     */
    explicit Matrix(std::string data);
    
    /**
     Constructor from a std::vector
     @param nRows Number of rows
     @param nCols Number of columns
     @param data Data vector
     */
    Matrix(std::size_t nRows, std::size_t nCols, std::vector<double> data);
    
    //////////////////////////////////////
    /// Accessors
    //////////////////////////////////////
    
    /**
     Number of rows of the matrix
     */
    std::size_t nRows() const
    {
        return m_nRows;
    }
    
    /**
     Number of columns of the matrix
     */
    std::size_t nCols() const
    {
        return m_nCols;
    }
    
    std::size_t Size() const
    {
        return m_nRows * m_nCols;
    }
    
    //////////////////////////////////////
    /// Operators
    //////////////////////////////////////
    
    /**
     Access a cell using parentheses
     @param i Row index
     @param j Column index
     */
    double& operator()(std::size_t i, std::size_t j)
    {
        // Using parentheses rather than square brackets is faster
        return m_data[offset(i, j)];
    }
    
    /**
     Access a cell using parentheses
     @param i Row index
     @param j Column index
     */
    const double& operator()(std::size_t i, std::size_t j) const
    {
        return m_data[offset(i, j)];
    }
    
    /**
     Addition
     @param rhs Matrix to add
     */
    Matrix& operator+=(const Matrix& rhs);
    
    /**
     Subtraction
     @param rhs Matrix to subtract
     */
    Matrix& operator-=(const Matrix& rhs);
    
    /**
     Scalar addition
     @param rhs Scalar to add to each element
     */
    Matrix& operator+=(double rhs);
    
    /**
     Scalar subtraction
     @param rhs Scalar to subtract to each element
     */
    Matrix& operator-=(double rhs);
    
    /**
     Scalar multiplication
     */
    Matrix& operator*=(double factor);
    
    /**
     Scalar division
     */
    Matrix& operator/=(double factor);
    
    ////////////////////////////////////////////
    /// Other useful methods
    ////////////////////////////////////////////
    
    /**
     Identity matrix
     @param dim Dimension
     */
    static Matrix Identity(std::size_t dim);
    
    /**
     Matrix transpose
     */
    Matrix transpose() const;
    
    /**
     Determinant of a square matrix
     */
    double Determinant() const;
    
    /**
     Inverse of a square matrix
     */
    Matrix Inverse() const;
    
    /**
     Asserts that two matrices are equal
     @param other Matrix to compare
     @param tolerance Elementwise tolerance level
     */
    void assertEquals(const Matrix& other, double tolerance) const;
    
    /**
     Extracts a given row from the matrix
     @param i Index of the row to extract
     */
    std::vector<double> ExtractRow(std::size_t i) const;
    
    /**
     Sets row of the matrix as input data vector
     @param i Index of the row to modify
     @param data Input data vector
     */
    void setRow(std::size_t i, const std::vector<double>& data);
    
    /**
     Converts a 1 x n Matrix into a row vector
     */
    std::vector<double> asRow() const;
    
    /**
     Extracts a given column from the matrix
     @param j Index of the column to extract
     */
    std::vector<double> ExtractColumn(std::size_t j) const;
    
    /**
     Sets column of the matrix as input data vector
     @param j Index of the column to modify
     @param data Input data vector
     */
    void setColumn(std::size_t j, const std::vector<double>& data);
    
    /**
     Converts a 1 x n Matrix into a column vector
     */
    std::vector<double> asColumn() const;
    
private:
    std::size_t m_nRows;
    std::size_t m_nCols;
    std::vector<double> m_data;
    
    /**
     Returns the vector position of the corresponding element
     @param i Row index
     @param j Column index
     */
    std::size_t offset(std::size_t i, std::size_t j) const
    {
        /*
         This allows us to store matrix in C-style, i.e. row-wise.
         Column-wise is Fortran-style.
         */
        ASSERT((i < m_nRows) && (j < m_nCols));
        return i * m_nCols + j;
    }
    
    /**
     Cofactor associated to the (i, j) element
     @param i Row index
     @param j Column index
     */
    Matrix Cofactor(std::size_t i, std::size_t j) const;
    
    /**
     Adjugate matrix
     */
    Matrix Adjugate() const;
};

/**
 Creates a matrix of zeros
 @param nRows Number of rows
 @param nCols Number of columns
 */
Matrix zeros(std::size_t nRows, std::size_t nCols);

/**
 Creates a matrix of ones
 @param nRows Number of rows
 @param nCols Number of columns
 */
Matrix ones(std::size_t nRows, std::size_t nCols);

/**
 Cholesky decomposition of input matrix
 @param m Input matrix
 */
Matrix Cholesky(const Matrix& m);

///////////////////////////////
///
/// Operators
///
///////////////////////////////

/**
 Operator adding two matrices of the same dimension
 @param lhs Left matrix
 @param rhs Right matrix
 */
Matrix operator+(const Matrix& lhs, const Matrix& rhs);

/**
 Operator adding a scalar to each element of the matrix
 @param lhs Matrix
 @param rhs Scalar to add
 */
Matrix operator+(const Matrix& lhs, double rhs);

/**
 Operator adding a matrix to a scalar
 @param lhs Scalar to add
 @param rhs Matrix
 */
inline Matrix operator+(double lhs, const Matrix& rhs)
{
    return rhs + lhs;
}

/**
 Operator substracting a matrix to another of the same dimension
 @param lhs Left matrix
 @param rhs Right matrix
 */
Matrix operator-(const Matrix& lhs, const Matrix& rhs);

/**
 Operator substracting a scalar to each element of the matrix
 @param lhs Matrix
 @param rhs Scalar to subtract
 */
Matrix operator-(const Matrix& lhs, double rhs);

/**
 Operator adding a scalar to the element-wise opposite of a matrix
 @param lhs Scalar
 @param rhs Matrix
 */
Matrix operator-(double lhs, const Matrix& rhs);

/**
 Matrix multiplication
 @param lhs Left Matrix
 @param rhs Right Matrix
 */
Matrix operator*(const Matrix& lhs, const Matrix& rhs);

/**
 Multiply a matrix by a scalar
 @param lhs Matrix to multiply
 @param rhs Multiplicative factor
 */
Matrix operator*(const Matrix& lhs, double rhs);

/**
 Multiply a matrix by a scalar
 @param lhs Multiplicative factor
 @param rhs Matrix to multiply
 */
inline Matrix operator*(double lhs, const Matrix& rhs)
{
    return rhs * lhs;
}

/**
 Matrix multiplication by a vector
 @param lhs Left Matrix
 @param rhs Right vector
 */
std::vector<double> operator*(const Matrix& lhs, const std::vector<double> rhs);

/**
 Matrix multiplication by a vector
 @param lhs Left Matrix
 @param rhs Right vector
 */
std::vector<double> operator*(const std::vector<double> lhs, const Matrix& rhs);

/**
 Prints the input matrix in the command line
 @param m Input matrix
 */
void printMatrix(const Matrix& m);

void testMatrix();

#endif /* Matrix_hpp */
