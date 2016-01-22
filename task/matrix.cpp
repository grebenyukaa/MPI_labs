#include <mpi.h>

#include <assert.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <omp.h>

#include "matrix.h"
#include "mpi_scope.h"

class Logger
{
public:
    Logger(const Matrix::index_type nodeid)
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

Matrix::Matrix(index_type cols, index_type rows)
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

/*const double Matrix::diagonality() const
{
    value_type norm = 0;
    for (index_type i = 0; i < m_rows; ++i)
        for (index_type j = i + 1; j < m_cols; ++j)
            norm += std::pow(at(i, j), 2);
    return norm;
}*/

void Matrix::compute_eigenvalues(const value_type& precision)
{
    value_type cur_norm;
    std::pair<index_type, index_type> ijmax = find_max_off_diagonal_norm(cur_norm);
    index_type imax = ijmax.first;
    index_type jmax = ijmax.second;
    value_type old_norm = cur_norm + 2 * precision;

    size_t iter;
    for (iter = 0; std::abs(cur_norm - old_norm) >= precision; ++iter)
    {
        //std::cout << *this << std::endl;
        jacoby_multiply(jmax, imax);
        //std::cout << *this << std::endl;

        old_norm = cur_norm;
        ijmax = find_max_off_diagonal_norm(cur_norm);
        imax = ijmax.first;
        jmax = ijmax.second;

        //if (iter % NTH_PRINT == 0)
        //    std::cout << "iteration " << iter << " delta = " << std::abs(cur_norm - old_norm) << std::endl;
        //std::cout << "--" << std::endl;
    }
    std::cout << "iteration " << iter << " delta = " << std::abs(cur_norm - old_norm) << std::endl;

    for (index_type i = 0; i < m_rows; ++i)
        m_eigenvalues[i] = at(i, i);
}

std::ostream& operator<<(std::ostream& o, const Matrix& m)
{
    for (Matrix::index_type i = 0; i < m.m_rows; ++i)
    {
        for (Matrix::index_type j = 0; j < m.m_cols; ++j)
        {
            o << std::setw(20) << m.at(i, j);
        }
        o << std::endl;
    }
    return o;
}

std::pair<Matrix::index_type, Matrix::index_type> Matrix::find_max_off_diagonal_norm_plain(value_type& norm)
{
    index_type imax = 0;
    index_type jmax = 0;
    value_type max = 0;
    for (index_type i = 0; i < m_rows; ++i)
        for (index_type j = i + 1; j < m_cols; ++j)
        {
            value_type absij = std::abs(at(i, j));
            norm += absij * absij;
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
Matrix::index_type find_abs_max_norm(const T* from, const size_t count, T& max_val, T& norm_part)
{
    using namespace std;

    norm_part = 0;
    max_val = 0;
    Matrix::index_type imax = 0;
    T abs_val = 0;

//#if defined(COMPUTATION_MPI_OMP)
//    #pragma omp parallel for private(abs_val)
//#endif
    for (size_t i = 0; i < count; ++i)
    {
        abs_val = std::abs(from[i]);
//#if defined(COMPUTATION_MPI_OMP)
//        #pragma omp atomic update
//#endif
        norm_part += abs_val * abs_val;
        if (max_val < abs_val)
        {
//#if defined(COMPUTATION_MPI_OMP)
//            #pragma omp atomic write
//#endif
            max_val = abs_val;
//#if defined(COMPUTATION_MPI_OMP)
//            #pragma omp atomic write
//#endif
            imax = i;
        }
    }

    return imax;
}

std::pair<Matrix::index_type, Matrix::index_type> Matrix::find_max_off_diagonal_norm_mpi_omp(value_type& norm)
{
    norm = 0;
    value_type max = 0;
    index_type imax = 0;
    index_type jmax = 0;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //int world_rank = MPI::COMM_WORLD.Get_rank();
    //int world_size = MPI::COMM_WORLD.Get_size();
    
    //Logger log(world_rank);
    try
    {
        assert(world_rank == 0 && "Should be executed in parent node!");
        //log "world_size = " << world_size << std::endl;

        index_type waiting_for = 0;
        for (index_type h = 0, wr = 0; h < m_rows; ++h, ++wr)
        {
            index_type row_sz = m_cols - h - 1;
            if (row_sz >= MIN_BATCH_SIZE)
            {
                index_type dest = wr % (world_size - 1) + 1;
                //log "sending row #" << h << " to " << dest << std::endl;
                {
                    MPI_Trace scp(MPI_Trace::ClSend);
                    MPI_Send(&m_data[h * m_rows + h + 1], m_cols - h - 1, MPI_DOUBLE, dest, h, MPI_COMM_WORLD);
                    //MPI::COMM_WORLD.Send(&m_data[h * m_rows + h + 1], m_cols - h - 1, MPI::DOUBLE, dest, h);
                }
                //log "  sent." << std::endl;
                ++waiting_for;
            }
            else
            {
                value_type norm_prt;
                value_type new_max;
                index_type new_jmax = find_abs_max_norm(&m_data[h * m_rows + h + 1], row_sz, new_max, norm_prt);
                norm += norm_prt;
                if (max < new_max)
                {
                    imax = h;
                    jmax = new_jmax + h + 1;
                    max = new_max;
                }
                //log "not sent row #" << h << ": size too small" << std::endl;
            }
        }

        //log "wating for reponses..." << std::endl;

        for (; waiting_for > 0; --waiting_for)
        {
            //log "[LEFT:" << std::setw(5) << waiting_for << "]" << " probing ..." << std::endl;

            MPI_Status status;
            //MPI::Status status;
            {
                volatile MPI_Trace scp(MPI_Trace::ClProbe); (void)scp;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                //MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);
            }

            index_type rowid = status.MPI_TAG/*status.Get_tag()*/;
            index_type nodeid = status.MPI_SOURCE/*status.Get_source()*/;
            //log "[LEFT:" << std::setw(5) << waiting_for << "] retrieving response from " << nodeid << " row #" << rowid << std::endl;

            std::vector<value_type> resp(3);
            {
                volatile MPI_Trace scp(MPI_Trace::ClRecv); (void)scp;
                MPI_Recv(&resp[0], 3, MPI_DOUBLE, nodeid, rowid, MPI_COMM_WORLD, &status);
                //MPI::COMM_WORLD.Recv(&resp[0], 3, MPI::DOUBLE, nodeid, rowid);
            }
            //log "  done." << std::endl;

            index_type new_jmax = (index_type)resp[0];
            value_type new_max = resp[1];
            value_type norm_prt = resp[2];

            norm += norm_prt;
            if (max < new_max)
            {
                imax = rowid;
                jmax = new_jmax + rowid + 1;
                max = new_max;
            }
        }

        //log "  done waiting." << std::endl;
    }
    catch (const std::exception& e)
    {
        //log e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, ABORTION_CODE);
        //MPI::COMM_WORLD.Abort(ABORTION_CODE);
    }

    return std::make_pair(imax, jmax);
}

void Matrix::find_max_mpi_server(const index_type cols, const index_type rows)
{
    (void)cols;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //int world_rank = MPI::COMM_WORLD.Get_rank();
    //int world_size = MPI::COMM_WORLD.Get_size();
    
    //Logger log(world_rank);
    try
    {
        assert(world_rank > 0 && "Should be executed in child node!");

        //index_type iter_cnt = (rows - MIN_BATCH_SIZE) % (world_size - 1) + 1;
        //log "node #" << world_rank << " is listening..." << std::endl;

        for (index_type iter = 0; /*iter < iter_cnt*/; ++iter)
        {
            MPI_Status status;
            //MPI::Status status;
            {
                volatile MPI_Trace scp(MPI_Trace::ClProbe); (void)scp;
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                //MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, status);
            }
            index_type rowid = status.MPI_TAG/*status.Get_tag()*/;
            //index_type nodeid = status.MPI_SOURCE/*status.Get_source()*/;

            int sz_to_read_;
            MPI_Get_count(&status, MPI_DOUBLE, &sz_to_read_);
            index_type sz_to_read = (index_type)sz_to_read_;
            //index_type sz_to_read = status.Get_count(MPI::DOUBLE);

            index_type row_sz = cols - rowid - 1;

            assert(sz_to_read == row_sz);

            std::vector<value_type> row(row_sz);
            {
                volatile MPI_Trace scp(MPI_Trace::ClRecv); (void)scp;
                MPI_Recv(&row[0], row_sz, MPI_DOUBLE, 0, rowid, MPI_COMM_WORLD, &status);
                //MPI::COMM_WORLD.Recv(&row[0], row_sz, MPI::DOUBLE, 0, rowid);
            }
            //log "[ITER: #" << std::setw(3) << iter << "] " << "received row #" << rowid << " from node #" << nodeid << " size = " << sz_to_read << std::endl;

            value_type new_max;
            value_type norm_prt;
            index_type new_jmax = find_abs_max_norm(&row[0], (index_type)row.size(), new_max, norm_prt);

            std::vector<value_type> resp(3);
            resp[0] = new_jmax;
            resp[1] = new_max;
            resp[2] = norm_prt;
            {
                volatile MPI_Trace scp(MPI_Trace::ClSend); (void)scp;
                MPI_Send(&resp[0], 3, MPI_DOUBLE, 0, rowid, MPI_COMM_WORLD);
                //MPI::COMM_WORLD.Send(&resp[0], 3, MPI::DOUBLE, 0, rowid);
            }
            //log "[ITER: #" << std::setw(3) << iter << "] " << "sent back answer (row #" << rowid << ")" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        //log e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, ABORTION_CODE);
        //MPI::COMM_WORLD.Abort(ABORTION_CODE);
    }
}

void Matrix::jacoby_multiply_plain(const index_type l, const index_type k)
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

#if defined(COMPUTATION_MPI_OMP)
    #pragma omp parallel for
#endif
    for (index_type h = 0; h < m_cols; ++h)
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
