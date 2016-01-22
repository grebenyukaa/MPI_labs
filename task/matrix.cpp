#include <assert.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <mpi.h>
#include <omp.h>

#include "matrix.h"

class Logger
{
public:
    Logger(const int nodeid)
    {
        std::ostringstream oss;
        oss << "node" << nodeid << ".log";
        m_ofs = new std::ofstream(oss.str().c_str(), std::ios_base::app);
    }

    ~Logger()
    {
        m_ofs->close();
        delete m_ofs;
    }

    template<class T>
    std::ofstream& operator<<(const T& data)
    {
        *m_ofs << data;
        return *m_ofs;
    }
private:
    std::ofstream* m_ofs;
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Matrix::Matrix(int cols, int rows)
    :
    m_cols(cols),
    m_rows(rows),
    m_data(rows * cols),
    m_eigenvalues(rows)
{}

Matrix::Matrix(const Matrix &other)
    :
    m_cols(other.m_cols),
    m_rows(other.m_rows),
    m_data(other.m_data)
{}

Matrix::~Matrix()
{}

const double Matrix::diagonality() const
{
    value_type norm = 0;
    for (int i = 0; i < m_rows; ++i)
        for (int j = i + 1; j < m_cols; ++j)
            norm += std::pow(at(i, j), 2);
    return norm;
}

void Matrix::compute_eigenvalues(const value_type& precision)
{
    value_type old_norm;
    value_type cur_norm = diagonality();
    int iter = 0;
    do
    {
        old_norm = cur_norm;

        std::pair<int, int> ijmax = find_max_off_diagonal();
        int imax = ijmax.first;
        int jmax = ijmax.second;

        //std::cout << *this << std::endl;
        jacoby_multiply(jmax, imax);
        //std::cout << *this << std::endl;

        cur_norm = diagonality();
        if (iter++ % NTH_PRINT == 0)
            std::cout << "iteration " << iter++ << " delta = " << std::abs(cur_norm - old_norm) << std::endl;
        //std::cout << "--" << std::endl;
    }
    while (std::abs(cur_norm - old_norm) >= precision);

    for (int i = 0; i < m_rows; ++i)
        m_eigenvalues[i] = at(i, i);
}

std::ostream& operator<<(std::ostream& o, const Matrix& m)
{
    for (int i = 0; i < m.m_rows; ++i)
    {
        for (int j = 0; j < m.m_cols; ++j)
        {
            o << std::setw(20) << m.at(i, j);
        }
        o << std::endl;
    }
    return o;
}

std::pair<int, int> Matrix::find_max_off_diagonal_plain()
{
    int imax = 0;
    int jmax = 0;
    value_type max = 0;
    for (int i = 0; i < m_rows; ++i)
        for (int j = i + 1; j < m_cols; ++j)
        {
            value_type absij = std::abs(at(i, j));
            if (max < absij)
            {
                max = absij;
                imax = i;
                jmax = j;
            }
        }
    return std::make_pair(imax, jmax);
}

template<class T>
int find_abs_max(const T* from, const int count, T& max_val)
{
    max_val = 0;
    int imax = 0;
#if defined(MATRIX_MUL_MPI_OMP)
    #pragma omp for
#endif
    for (int i = 0; i < count; ++i)
    {
        T abs_val = std::abs(from[i]);
        if (max_val < abs_val)
        {
#if defined(MATRIX_MUL_MPI_OMP)
            #pragma omp critical (max_update)
#endif
            max_val = abs_val;
            imax = i;
        }
    }
    return imax;
}

std::pair<int, int> Matrix::find_max_off_diagonal_mpi_omp()
{
    value_type max = 0;
    int imax = 0;
    int jmax = 0;

    int world_rank = MPI::COMM_WORLD.Get_rank();
    int world_size = MPI::COMM_WORLD.Get_size();
    //Logger log(world_rank);
    try
    {
        assert(world_rank == 0 && "Should be executed in parent node!");
        //log << "world_size = " << world_size << std::endl;

        int waiting_for = 0;
        for (int h = 0, wr = 0; h < m_rows; ++h, ++wr)
        {
            int row_sz = m_cols - h - 1;
            if (row_sz >= MIN_BATCH_SIZE)
            {
                int dest = wr % (world_size - 1) + 1;
                //log << "sending row #" << h << " to " << dest << std::endl;
                MPI::COMM_WORLD.Send(&m_data[h * m_rows + h + 1], m_cols - h - 1, MPI::DOUBLE, dest, h);
                //log << "  sent." << std::endl;
                ++waiting_for;
            }
            else
            {
                value_type new_max;
                int new_jmax = find_abs_max(&m_data[h * m_rows + h + 1], row_sz, new_max);
                if (max < new_max)
                {
                    imax = h;
                    jmax = new_jmax + h + 1;
                    max = new_max;
                }
                //log << "not sent row #" << h << ": size too small" << std::endl;
            }
        }

        //log << "wating for reponses..." << std::endl;

        for (; waiting_for > 0; --waiting_for)
        {
            //log << "[LEFT:" << std::setw(5) << waiting_for << "]" << " probing ..." << std::endl;

            MPI::Status status;
            MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);

            int rowid = status.Get_tag();
            int nodeid = status.Get_source();
            //log << "[LEFT:" << std::setw(5) << waiting_for << "] retrieving response from " << nodeid << " row #" << rowid << std::endl;

            std::vector<value_type> resp(2);
            MPI::COMM_WORLD.Recv(&resp[0], 2, MPI::DOUBLE, nodeid, rowid);
            //log << "  done." << std::endl;

            int new_jmax = (int)resp[0];
            value_type new_max = resp[1];
            if (max < new_max)
            {
                imax = rowid;
                jmax = new_jmax + rowid + 1;
                max = new_max;
            }
        }

        //log << "  done waiting." << std::endl;
    }
    catch (const std::exception& e)
    {
        //log << e.what() << std::endl;
        MPI::COMM_WORLD.Abort(ABORTION_CODE);
    }

    return std::make_pair(imax, jmax);
}

void Matrix::find_max_mpi_server(const int cols, const int rows)
{
    (void)cols;

    int world_rank = MPI::COMM_WORLD.Get_rank();
    int world_size = MPI::COMM_WORLD.Get_size();
    //Logger log(world_rank);
    try
    {
        assert(world_rank > 0 && "Should be executed in child node!");

        int iter_cnt = (rows - MIN_BATCH_SIZE) % (world_size - 1) + 1;
        //log << "node #" << world_rank << " is going to process up to " << iter_cnt << " requests..." << std::endl;

        for (int iter = 0; /*iter < iter_cnt*/; ++iter)
        {
            MPI::Status status;
            MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, status);
            int rowid = status.Get_tag();
            int nodeid = status.Get_source();
            int sz_to_read = status.Get_count(MPI::DOUBLE);
            int row_sz = cols - rowid - 1;

            assert(sz_to_read == row_sz);

            std::vector<value_type> row(row_sz);
            MPI::COMM_WORLD.Recv(&row[0], row_sz, MPI::DOUBLE, 0, rowid);
            //log << "[ITER: #" << std::setw(3) << iter << "] " << "received row #" << rowid << " from node #" << nodeid << " size = " << sz_to_read << std::endl;

            value_type new_max;
            int new_jmax = find_abs_max(&row[0], (int)row.size(), new_max);

            std::vector<value_type> resp(2);
            resp[0] = new_jmax;
            resp[1] = new_max;
            MPI::COMM_WORLD.Send(&resp[0], 2, MPI::DOUBLE, 0, rowid);
            //log << "[ITER: #" << std::setw(3) << iter << "] " << "sent back answer (row #" << rowid << ")" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        //log << e.what() << std::endl;
        MPI::COMM_WORLD.Abort(ABORTION_CODE);
    }
}

void Matrix::jacoby_multiply_plain(const int l, const int k)
{
    using namespace std;

    assert(m_cols == m_rows);
    if (abs(at(k, l)) < 1e-12) return;

    value_type beta = (at(l, l) - at(k, k)) / at(k, l) / 2;
    value_type t = sgn(beta) / (abs(beta) + sqrt(beta*beta + 1));
    value_type c = 1 / sqrt(t*t + 1);
    value_type s = c*t;
    value_type ro = s / (1 + c);

    value_type kk = at(k, k) - t * at(k, l);
    value_type ll = at(l, l) + t * at(k, l);

#if defined(MATRIX_MUL_MPI_OMP)
    #pragma omp for
#endif
    for (int h = 0; h < m_cols; ++h)
    {
        if (h == k || h == l) continue;
        value_type hk = at(h, k) - s*(at(h, l) + ro*at(h, k));
        value_type hl = at(h, l) + s*(at(h, k) - ro*at(h, l));
        at(h, k) = hk;
        at(k, h) = hk;
        at(h, l) = hl;
        at(l, h) = hl;
    }
    at(k, k) = kk;
    at(l, l) = ll;
    at(k, l) = 0;
    at(l, k) = 0;
}
