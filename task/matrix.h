#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

#define NTH_PRINT 1

#define MIN_BATCH_SIZE 50
#define ABORTION_CODE -1

//#define COMPUTATION_PLAIN
//#define COMPUTATION_MPI
//#define COMPUTATION_MPI_OMP

class Matrix
{
public:
    typedef double value_type;
    typedef unsigned long long index_type;

    Matrix(index_type cols, index_type rows);
    //rotation matrix
    Matrix(index_type cols, index_type rows, index_type col, index_type row, value_type s, value_type c);
    Matrix(const Matrix& other);
    virtual ~Matrix();

    void compute_eigenvalues(const value_type& precision);
    //const double diagonality() const;

    inline value_type& at(const index_type& i, const index_type& j) { return m_data[i * m_rows + j]; }
    inline const value_type& at(const index_type& i, const index_type& j) const { return m_data[i * m_rows + j]; }

    inline const index_type getColCount() const { return m_cols; }
    inline const index_type getRowCount() const { return m_rows; }
    inline const std::vector<value_type>& getEigenValues() const { return m_eigenvalues; }

    friend std::ostream& operator<<(std::ostream& o, const Matrix& m);

    static void find_max_mpi_server(const index_type cols, const index_type rows);
private:
    //seemed to be a good idea, but really nothing to parallelize here
    inline void jacoby_multiply(const index_type col, const index_type row) { jacoby_multiply_plain(col, row); }
    void jacoby_multiply_plain(const index_type col, const index_type row);

    //parallelizing max counting
    inline std::pair<index_type, index_type> find_max_off_diagonal_norm(value_type& norm)
    {
        #if   defined(COMPUTATION_PLAIN)
            return find_max_off_diagonal_norm_plain(norm);
        #elif defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
            return find_max_off_diagonal_norm_mpi_omp(norm);
        #endif
    }

    std::pair<index_type, index_type> find_max_off_diagonal_norm_plain(value_type& norm);
    std::pair<index_type, index_type> find_max_off_diagonal_norm_mpi(value_type& norm);
    std::pair<index_type, index_type> find_max_off_diagonal_norm_mpi_omp(value_type& norm);

    index_type m_cols;
    index_type m_rows;
    std::vector<value_type> m_data;
    std::vector<value_type> m_eigenvalues;
};

#endif //MATRIX_H
