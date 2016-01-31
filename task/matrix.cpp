#include <mpi.h>

#include <assert.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <omp.h>

#include "matrix.h"
#include "mpi_scope.h"
#include "utils.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//

Matrix::Matrix(index_type size)
    :
    m_cols(size),
    m_rows(size),
    m_data(APUtils::sum(1, size)),
    m_eigenvalues(size)
{}

Matrix::Matrix(const Matrix &other)
    :
    m_cols(other.m_cols),
    m_rows(other.m_rows),
    m_data(other.m_data)
{}

Matrix::~Matrix()
{}

Matrix::value_type& Matrix::at(index_type i, index_type j)
{
    assert(j != i);
    if (j < i) std::swap(i, j);
    index_type sum = APUtils::sum(m_rows - 1 - i, m_rows - 1);
    return m_data[sum - (m_rows - 1) + j - 1];
}

const Matrix::value_type& Matrix::at(index_type i, index_type j) const
{
    assert(j != i);
    if (j < i) std::swap(i, j);
    index_type sum = APUtils::sum(m_rows - 1 - i, m_rows - 1);
    return m_data[sum - (m_rows - 1) + j - 1];
}

Matrix::value_type& Matrix::at_diag(index_type i)
{
    return m_eigenvalues[i];
}

const Matrix::value_type& Matrix::at_diag(index_type i) const
{
    return m_eigenvalues[i];
}

std::ostream& operator<<(std::ostream& o, const Matrix& m)
{
    for (Matrix::index_type i = 0; i < m.m_rows; ++i)
    {
        o << std::setw(20 * (i + 1)) << m.m_eigenvalues[i];
        for (Matrix::index_type j = i + 1; j < m.m_cols; ++j)
        {
            o << std::setw(20) << m.at(i, j);
        }
        o << std::endl;
    }
    return o;
}

void Matrix::compute_eigenvalues(const value_type& precision)
{
    value_type cur_norm;
    std::pair<index_type, index_type> ijmax = find_max_off_diagonal_norm(cur_norm);
    index_type imax = ijmax.first;
    index_type jmax = ijmax.second;
    value_type old_norm = 0;
    //std::cout << "iteration " << -1 << " delta = " << std::abs(cur_norm - old_norm) << " max at (" << imax << ", " << jmax << ") = " << std::setw(20) << at(imax, jmax) << std::endl;

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
        //    std::cout << "iteration " << iter << " delta = " << std::abs(cur_norm - old_norm) << " max at (" << imax << ", " << jmax << ") = " << std::setw(20) << at(imax, jmax) << std::endl;
        //std::cout << "--" << std::endl;

#ifdef MPITV_ABORT_AFTER
        if (iter == MPITV_ABORT_AFTER)
        {
            std::cout << "mpitv: aborting with code " << MPITV_ABORTION_CODE << std::endl;
            //MPI_Finalize();
            //exit(MPITV_ABORTION_CODE);
            MPI_Pcontrol(MPI_TRACEFLUSH);
            MPI_Abort(MPI_COMM_WORLD, MPITV_ABORTION_CODE);
            return;
        }
#endif
    }
    std::cout << "iteration " << iter << " delta = " << std::abs(cur_norm - old_norm) << std::endl;
}

//

std::pair<Matrix::index_type, Matrix::index_type> Matrix::flat_idx_to_pair(index_type idx)
{
    int i = 0;
    int j = 0;
    int rowsize = m_cols - 1;
    for (; (idx >= rowsize); idx -= rowsize--)
        ++i;
    assert(idx >= 0);
    j = i + idx + 1;
    return std::make_pair(i, j);
}

// plain

std::pair<Matrix::index_type, Matrix::index_type> Matrix::find_max_off_diagonal_norm_plain(value_type& norm)
{
    index_type imax = 0;
    value_type max = 0;
    for (size_t i = 0; i < m_data.size(); ++i)
    {
        value_type absv = std::abs(m_data[i]);
        norm += absv * absv;
        if (max < absv)
        {
            max = absv;
            imax = i;
        }
    }
    return flat_idx_to_pair(imax);
}

void Matrix::jacoby_multiply_plain(const index_type l, const index_type k)
{
    using namespace std;

    assert(k != l);
    value_type kl = at(k, l);
    if (abs(kl) < 1e-12) return;

    value_type beta = (at_diag(l) - at_diag(k)) / kl / 2;
    value_type t = sgn(beta) / (abs(beta) + sqrt(beta*beta + 1));
    value_type c = 1 / sqrt(t*t + 1);
    value_type s = c*t;
    value_type ro = s / (1 + c);

    at_diag(k) -= t * kl;
    at_diag(l) += t * kl;

#if defined(COMPUTATION_MPI_OMP)
    #pragma omp parallel for
#endif
    for (index_type h = 0; h < m_rows; ++h)
    {
        if (h == l || h == k) continue;
        //std::cout << "(" << h << ", " << k << ")" << std::endl;
        //std::cout << "(" << h << ", " << l << ")" << std::endl;
        value_type hk = at(h, k) - s*(at(h, l) + ro*at(h, k));
        value_type hl = at(h, l) + s*(at(h, k) - ro*at(h, l));
        //at(h, k) = hk;
        at(k, h) = hk;
        //at(h, l) = hl;
        at(l, h) = hl;
    }
    at(k, l) = 0;
    //at(l, k) = 0;
}

//mpi, mpi_omp

template<class T>
Matrix::index_type find_abs_max_norm(const T* from, const size_t count, T& max_val, T& norm_part)
{
    using namespace std;

    norm_part = 0;
    max_val = 0;
    Matrix::index_type imax = 0;
    T abs_val = 0;
    size_t i;

#if defined(COMPUTATION_MPI_OMP)
//    #pragma omp parallel for private(abs_val, i)
#endif
    for (i = 0; i < count; ++i)
    {
        abs_val = std::abs(from[i]);
#if defined(COMPUTATION_MPI_OMP)
//        #pragma omp atomic update
#endif
        norm_part += abs_val * abs_val;
        if (max_val < abs_val)
        {
#if defined(COMPUTATION_MPI_OMP)
//            #pragma omp atomic write
#endif
            max_val = abs_val;
#if defined(COMPUTATION_MPI_OMP)
//            #pragma omp atomic write
#endif
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

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    Logger log(world_rank);
    try
    {
        assert(world_rank == 0 && "Should be executed in parent node!");
        log << "world_size = " << world_size << std::endl;

        size_t batch_sz = m_data.size() % (world_size - 1) == 0 ? m_data.size() / (world_size - 1) : m_data.size() / (world_size - 1) + 1;
        for (int node_id = 1; node_id < world_size; ++node_id)
        {
            {
                Matrix::index_type offset = (node_id - 1) * batch_sz;
                //volatile MPI_Trace scp(MPI_Trace::ClSend);
                MPI_TRACE_EVENT(MPI_Send(&m_data[offset], std::min(batch_sz, m_data.size() - offset), MPI_DOUBLE, node_id, node_id, MPI_COMM_WORLD), MPI_Trace::ClSend);
            }
        }

        for (int node_id = 1; node_id < world_size; ++node_id)
        {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            index_type nodeid = status.MPI_SOURCE;
            index_type tag = status.MPI_TAG;

            std::vector<value_type> resp(3);
            {
                //volatile MPI_Trace scp(MPI_Trace::ClRecv);
                MPI_TRACE_EVENT(MPI_Recv(&resp[0], 3, MPI_DOUBLE, nodeid, tag, MPI_COMM_WORLD, &status), MPI_Trace::ClRecv);
            }

            index_type new_imax = (nodeid - 1) * batch_sz + (index_type)resp[0];
            value_type new_max = resp[1];
            value_type norm_prt = resp[2];

            norm += norm_prt;
            if (max < new_max)
            {
                imax = new_imax;
                max = new_max;
            }
        }
    }
    catch (const std::exception& e)
    {
        //log << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, ABORTION_CODE);
    }

    return flat_idx_to_pair(imax);
}

void Matrix::find_max_mpi_server(const index_type cols, const index_type rows)
{
    (void)cols;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    Logger log(world_rank);
    try
    {
        assert(world_rank > 0 && "Should be executed in child node!");

        //log << "node #" << world_rank << " is listening..." << std::endl;

        for (index_type iter = 0; /*iter < iter_cnt*/; ++iter)
        {
            MPI_Status status;
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            index_type tag = status.MPI_TAG;

            index_type sz_to_read;
            MPI_Get_count(&status, MPI_DOUBLE, &sz_to_read);

            std::vector<value_type> data(sz_to_read);
            {
                //volatile MPI_Trace scp(MPI_Trace::ClRecv);
                MPI_TRACE_EVENT(MPI_Recv(&data[0], sz_to_read, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status), MPI_Trace::ClRecv);
            }

            value_type new_max;
            value_type norm_prt;
            index_type new_jmax = find_abs_max_norm(&data[0], (index_type)data.size(), new_max, norm_prt);

            std::vector<value_type> resp(3);
            resp[0] = new_jmax;
            resp[1] = new_max;
            resp[2] = norm_prt;
            {
                //volatile MPI_Trace scp(MPI_Trace::ClSend);
                MPI_TRACE_EVENT(MPI_Send(&resp[0], 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD), MPI_Trace::ClSend);
            }
        }
    }
    catch (const std::exception& e)
    {
        //log << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, ABORTION_CODE);
    }
}
