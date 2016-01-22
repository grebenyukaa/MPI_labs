#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

#define MIN_BATCH_SIZE 5
#define ABORTION_CODE -1

//#define MATRIX_MUL_PLAIN
#define MATRIX_MUL_MPI

class Matrix
{
public:
    typedef double value_type;

    Matrix(int cols, int rows);
    //rotation matrix
    Matrix(int cols, int rows, int col, int row, value_type s, value_type c);
    Matrix(const Matrix& other);
    virtual ~Matrix();

    void compute_eigenvalues(const value_type& precision);
    const double diagonality() const;

    inline value_type& at(const int& i, const int& j) { return m_data[i * m_rows + j]; }
    inline const value_type& at(const int& i, const int& j) const { return m_data[i * m_rows + j]; }

    inline const int getColCount() const { return m_cols; }
    inline const int getRowCount() const { return m_rows; }
    inline const std::vector<value_type>& getEigenValues() const { return m_eigenvalues; }

    friend std::ostream& operator<<(std::ostream& o, const Matrix& m);

    static void find_max_mpi_server(const int cols, const int rows);
private:
    //seemed to be a good idea, but really nothing to parallelize here
    inline void jacoby_multiply(const int col, const int row) { jacoby_multiply_plain(col, row); }

    //parallelizing max counting
    inline std::pair<int, int> find_max_off_diagonal()
    {
        #if   defined(MATRIX_MUL_PLAIN)
            return find_max_off_diagonal_plain();
        #elif defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
            return find_max_off_diagonal_mpi_omp();
        #endif
    }

    void jacoby_multiply_plain(const int col, const int row);

    std::pair<int, int> find_max_off_diagonal_plain();
    std::pair<int, int> find_max_off_diagonal_mpi();
    std::pair<int, int> find_max_off_diagonal_mpi_omp();

    int m_cols;
    int m_rows;
    std::vector<value_type> m_data;
    std::vector<value_type> m_eigenvalues;
};

#endif //MATRIX_H
