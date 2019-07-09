#include <hpx/hpx_init.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/include/compute.hpp>
#include <hpx/include/parallel_transform.hpp>

#include <boost/program_options.hpp>

#include <algorithm>
#include <vector>
#include <iostream>

#define CALC_TYPE double

std::vector<CALC_TYPE> read_(std::string filename, std::size_t row, std::size_t col, std::size_t N, std::size_t T)
{
   std::ostringstream os;
   os << "read r:" << row << " c:" << col << " N:" << N << std::endl;
   std::cout << os.str();

   std::vector<CALC_TYPE> v;
   v.resize(N * N);
   std::srand(row * T + col);
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
         v[i * N + j] = (double) (std::rand() % 10000 - 5000) / 1000;
         std::ostringstream os;
         os << filename << " " << row * N + i << " " << col * N + j << " " << v[i * N + j] << std::endl;
         std::cout << os.str();
      }
   }
   return v;
}

//return inv
std::vector<CALC_TYPE> inversion(hpx::shared_future<std::vector<CALC_TYPE>> ft_A,
                           std::size_t k, std::size_t N)
{
   std::ostringstream os;
   os << "inv k:" << k << " N:" << N << std::endl;
   std::cout << os.str();
   auto A = ft_A.get();

   double tmp;
   std::vector<CALC_TYPE> v;
   v.resize(N * N);

   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
         v[i * N + j] = 0;
      }
      v[i * N + i] = 1;
   }
   for(int k = 0; k < N; ++k) {
      tmp = A[k * N + k];
      for(int j = 0; j < N; ++j) {
         v[k * N + j] /= tmp;
         A[k * N + j] /= tmp;
      }

      // can be parallalized
      for(int i = 0; i < k; ++i) {
         tmp = A[i * N + k];
         for(int j = 0; j < N; ++j) {
            v[i * N + j] -= tmp * v[k * N + j];
            A[i * N + j] -= tmp * A[k * N + j];
         }
      }
      for(int i = k + 1; i < N; ++i) {
         tmp = A[i * N + k];
         for(int j = 0; j < N; ++j) {
            v[i * N + j] -= tmp * v[k * N + j];
            A[i * N + j] -= tmp * A[k * N + j];
         }
      }
   }
   return v;
}

//A = A * B
std::vector<CALC_TYPE> pmm(hpx::shared_future<std::vector<CALC_TYPE>> ft_A,
                           hpx::shared_future<std::vector<CALC_TYPE>> ft_B,
                           std::size_t row, std::size_t col, std::size_t N)
{
   std::ostringstream os;
   os << "pmm r:" << row << " c:" << col << " N:" << N << std::endl;
   std::cout << os.str();

   std::vector<CALC_TYPE> v;
   v.resize(N * N);
   auto A = ft_A.get();
   auto B = ft_B.get();

   // i and j can be parallalized
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
         v[i * N + j] = 0;
         for(int k = 0; k < N; ++k) {
            v[i * N + j] += A[i * N + k] * B[k * N + j];
         }
      }
   }
   return v;
}

//C = C - A * B
std::vector<CALC_TYPE> pmm_d(hpx::shared_future<std::vector<CALC_TYPE>> ft_A,
                             hpx::shared_future<std::vector<CALC_TYPE>> ft_B,
                             hpx::shared_future<std::vector<CALC_TYPE>> ft_C,
                             std::size_t step, std::size_t row, std::size_t col, std::size_t N)
{
   std::ostringstream os;
   os << "pmm_d step" << step << " r:" << row << " c:" << col << " N:" << N << std::endl;
   std::cout << os.str();

   auto A = ft_A.get();
   auto B = ft_B.get();
   auto C = ft_C.get();

   // i and j can be parallalized
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
         for(int k = 0; k < N; ++k) {
            C[i * N + j] -= A[i * N + k] * B[k * N + j];
         }
      }
   }
   return C;
}

void write_(std::vector<CALC_TYPE> A, std::string filename, std::size_t row, std::size_t col, std::size_t N)
{
   std::ostringstream os;
   os << "write r:" << row << " c:" << col << " N:" << N << std::endl;
   std::cout << os.str();
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
         std::ostringstream os;
         os << filename << " " << row * N + i << " " << col * N + j << " " << A[i * N + j] << std::endl;
         std::cout << os.str();
      }
   }
}

void lu_tiled(std::size_t N, std::size_t T)
{
    std::vector<hpx::shared_future<std::vector<CALC_TYPE>>> ft_tiles;
    ft_tiles.resize(T * T);
    std::vector<hpx::shared_future<std::vector<CALC_TYPE>>> ft_inv;
    ft_inv.resize(T);

    for (std::size_t i = 0; i < T; ++i)
    {
       for (std::size_t j = 0; j < T; ++j)
       {
          ft_tiles[i * T + j] = hpx::async(&read_, "a.bin", i, j, N, T);
       }
    }

    for (std::size_t k = 0; k < T - 1; ++k)
    {
       ft_inv[k] = hpx::async(&inversion, ft_tiles[k * T + k], k, N);
       for (std::size_t i = k + 1; i < T; ++i)
       {
          ft_tiles[i * T + k] = hpx::async(&pmm, ft_tiles[i * T + k], ft_inv[k], i, k, N);
          for (std::size_t j = k + 1; j < T; ++j)
          {
             ft_tiles[i * T + j] = hpx::async(&pmm_d, ft_tiles[i * T + k], ft_tiles[k * T + j], ft_tiles[i * T + j], k, i, j, N);
          }
       }
    }

    for (std::size_t i = 0; i < T; ++i)
    {
       for (std::size_t j = 0; j < T; ++j)
       {
          write_(ft_tiles[i * T + j].get(), "lu.bin", i, j, N);
       }
    }
}

int hpx_main(boost::program_options::variables_map& vm)
{
    std::size_t N = vm["N"].as<std::size_t>();
    std::size_t T = vm["T"].as<std::size_t>();

    hpx::util::high_resolution_timer t;
    lu_tiled(N, T);
    double elapsed = t.elapsed();
    std::cout << "Elapsed " << elapsed << " s\n";
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    using namespace boost::program_options;

    options_description desc_commandline;
    desc_commandline.add_options()
        ("N", value<std::size_t>()->default_value(10),
         "Dimension of each Tile (N*N elements per tile)")
        ("T", value<std::size_t>()->default_value(10),
         "Number of Tiles in each dimension (T*T tiles)")
    ;

    return hpx::init(desc_commandline, argc, argv);
}
