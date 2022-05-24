#ifndef PARALLEL_MEDIAN_MEDIAN_HPP
#define PARALLEL_MEDIAN_MEDIAN_HPP
#include <vector>
#include <mpi.h>
#include <random>
#include <iostream>
#include <iterator>
MPI_Comm   new_comm;
using namespace std;
template<class T> void split_vec(std::vector<T>& vec, T pivot, std::vector<T>& le, std::vector<T>& gr, std::vector<T>& pi){

    auto m_it1 = std::partition(std::begin(vec), std::end(vec), [pivot](const T& em){ return em < pivot; });
    auto m_it2 = std::partition(m_it1, std::end(vec), [pivot](const T& em){ return em == pivot; });
    //auto m_it3 = std::partition(m_it2, std::end(vec), [pivot](const T& em){ return !(em > pivot); });

    le.assign(std::begin(vec), m_it1);
    pi.assign(m_it1, m_it2);
    gr.assign(m_it2, std::end(vec));
}

template<class T>
std::tuple<typename std::vector<T>::iterator, typename std::vector<T>::iterator>
split_vec(std::vector<T>& vec, T pivot){
    auto m_it1 = std::partition(std::begin(vec), std::end(vec), [pivot](const T& em){ return em < pivot; });
    auto m_it2 = std::partition(m_it1, std::end(vec), [pivot](const T& em){ return em == pivot; });
    //auto m_it3 = std::partition(m_it2, std::end(vec), [pivot](const auto& em){ return !(em > pivot); });
    return std::make_tuple(m_it1,m_it2);
}

template<class T> double nlogn_median(std::vector<T> v){
    std::sort(v.begin(), v.end());
    if(v.size() % 2)
        return v.at((v.size() / 2.0));
    else
        return 0.5 * (v.at(v.size() / 2 - 1) + v.at(v.size() / 2));
}
template<class T> double median(std::vector<T> x, size_t look_for) {
    std::random_device rd;
    std::vector<T> le, gr, pi;
    do {
        T pivot = x.at(rd() % x.size());
        split_vec(x, pivot, le, gr, pi);
        if(look_for < le.size()) {
            x = le;
        } else if (look_for < le.size() + pi.size()) {
            return pivot;
        } else {
            x = gr;
            look_for = look_for - le.size() - pi.size();
        }
    } while(true);
}
template<class T> double median(std::vector<T> x){
    if(x.size() % 2){
        return median(x, x.size() / 2);
    } else {
        return 0.5 * (median(x, x.size() / 2 - 1) + median(x, x.size() / 2));
    }
}
namespace par {
    int ws,rk;
    template<class T>
    MPI_Datatype get_mpi_type(){
        return MPI_UNSIGNED_LONG_LONG;
    }
    namespace {
        template<class T>
        double median(std::vector<T> x, size_t look_for) {
            std::random_device rd;
            if (MPI_COMM_NULL != new_comm) {
                MPI_Comm_size(new_comm, &ws);
                MPI_Comm_rank(new_comm, &rk);
            }
            /**
             * Until vec.size() == 1
             * 1. Someone selects pivot and broadcast
             * 2. All: split into le and gr
             * 3. All: keeps le if SUM(le_p) > SUM(gr_p) for all p in procs
             * 4. Go back to 1.
             */
            std::array<size_t, 3> split_sizes{};
            size_t size, total_size, lb, ub, ipivot;
            T pivot;
            do {
                size = x.size();

                MPI_Reduce(&size, &total_size, 1, get_mpi_type<size_t>(), MPI_SUM, 0, new_comm);
                MPI_Scan(&size, &ub, 1, get_mpi_type<size_t>(), MPI_SUM, new_comm);
                lb = ub - size;
                ipivot = (rd() % total_size);

                MPI_Bcast(&ipivot, 1, get_mpi_type<size_t>(), 0, new_comm);
                if(lb <= ipivot && ipivot < ub) {
                    pivot = x.at(ipivot - lb);
                    for(auto pe = 0; pe < ws; ++pe)
                        if(pe != rk) MPI_Send(&pivot, 1, get_mpi_type<T>(), pe, 999, new_comm);
                } else {
                    MPI_Recv(&pivot, 1, get_mpi_type<T>(), MPI_ANY_SOURCE, 999, new_comm, MPI_STATUSES_IGNORE);
                }
                auto t = split_vec(x, pivot);
                auto lit = std::get<0>(t);
                auto pit = std::get<1>(t);

                size_t le_size = std::distance(std::begin(x), lit),
                        pi_size = std::distance(lit, pit);

                split_sizes = {le_size, pi_size};

                MPI_Allreduce(MPI_IN_PLACE, &split_sizes, 2, get_mpi_type<size_t>(), MPI_SUM, new_comm);

                if(look_for < split_sizes[0]) {
                    x = std::vector<T>(std::begin(x), lit);
                } else if (look_for < split_sizes[0] + split_sizes[1]) {
                    return pivot;
                } else {
                    x = std::vector<T>(pit, std::end(x));
                    look_for = look_for - split_sizes[0] - split_sizes[1];
                }
            } while(true);
        }
    }
    template<class T>
    double median(const std::vector<T>& x) {
        size_t total_size, size = x.size();
        MPI_Allreduce(&size, &total_size, 1, get_mpi_type<size_t>(), MPI_SUM, new_comm);
        if(total_size % 2) {
            return median<T>(x, total_size / 2);
        } else {
            return 0.5 * (median<T>(x, total_size / 2 - 1) + median<T>(x, total_size / 2));
        }
    }
    template<class InputIterator>
    double median(InputIterator beg, InputIterator end) {
        std::vector<typename std::iterator_traits<InputIterator>::value_type> x(beg, end);
        return median(x);
    }
}
#endif //PARALLEL_MEDIAN_MEDIAN_HPP