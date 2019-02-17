/* copyright Grant Rostig (c)2019 see LICENSE file */
#include <iostream>
#include <iomanip>
#include <string>
#include <string_view>

#include <map>
#include <unordered_map>
#include <iterator>
#include <random>
#include <thread>
#include <algorithm>
#include <functional>
#include <future>
#include <cassert>
#include <stdexcept>
#include <boost/math/special_functions/binomial.hpp>
//#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <tr1/cmath>  // to pick up beta()
#include <eigen3/Eigen/Dense>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;

namespace grostig  {

//   * Bjarne Stroustrup. The C++ Programming Language (4th edition). 2013.
//   ISBN: 978-0-321-56384-2. Chapter 11.5: Explicit type conversion. page 299.
template <class Target, class Source>
Target narrow_cast_runtime(Source v)
{
  auto r = static_cast<Target>(v);  // todo: probably some undefined behavior with float/double
  if (static_cast<Source>(r)!=v)
    throw std::runtime_error("narrow_cast<>() failed");
  return r;
}

class thread_guard
{
    std::thread& t;
public:
    explicit thread_guard(std::thread& t_):
        t(t_)
    {}
    ~thread_guard()
    {
        if(t.joinable())
        {
            t.join();
        }
    }
    thread_guard(thread_guard const&)=delete;
    thread_guard& operator=(thread_guard const&)=delete;
};

/* std::cout << "Pascal's triangle:\n";  // may be wrong.
for(int n = 1; n < 10; ++n) {
    std::cout << std::string(20-n*2, ' ');
    for(int k = 1; k < n; ++k)
        std::cout << std::setw(3) << binomial_coefficient_nCk(n,k) << ' ';
    std::cout << '\n';
}
double binomial_coefficient_nCk(int const n, int const k) {  // THIS IS partly WRONG, use boost.
    if (n == 0 || n == k) return 1;
    return 1/ ((n+1) * std::tr1::beta( n-k+1, k+1 ));
} */

}
std::string      PGM_NAME               {"binomial"};
constexpr int    NUM_CPU_CORES =        8;
constexpr int    MAX_THREADS =          NUM_CPU_CORES;
constexpr int    MAX_ASYNC_FUTURES =    NUM_CPU_CORES;
constexpr size_t MAX_BINOMIAL_TRIALS =  1000;
constexpr size_t TRIALS_ARRAY_SZ =      (MAX_BINOMIAL_TRIALS+1)*MAX_THREADS;
using Histogram = std::array< size_t, MAX_BINOMIAL_TRIALS >;     // todo : gets quite large...
using Histogram_Parallel = std::array< size_t, TRIALS_ARRAY_SZ >;     // todo : gets quite large...

void boost_log_init() {  // BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
        boost::log::core::get()->set_filter( boost::log::trivial::severity >= boost::log::trivial::info);
}

double binomial_distribution_probability_mass_function( size_t const n_trials, size_t const k_successes, double const probability_of_success) {
    double bc =           boost::math::binomial_coefficient<double>( 5, 2);
    double p_k =          std::pow(probability_of_success, k_successes);
    errno = 0;
    double p_complement = std::pow(1 - probability_of_success, n_trials - k_successes);
    if ( 0 != errno ) { std::perror("binomial:ERROR: pow() failed."); }  // todo: what type do I need to use PGM_NAME
    return bc * p_k * p_complement;         // todo:  what happens if this cal overflows?  How do we check for that?
}

double print_binomial_dist_PMF(size_t const n_samples_drawn, size_t const k_successes, double const mu_probability_in_bin)
{
    try {
    auto r1 = binomial_distribution_probability_mass_function(n_samples_drawn, k_successes,  mu_probability_in_bin);
    cout << "binomial_probability_mass_function(n_trials,k_successes,p): " << n_samples_drawn << ", " << k_successes << ", " << mu_probability_in_bin << endl;
    cout << r1 << endl;
    return r1;
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

void binomial_rel_freq_histogram_map(size_t const trials, double const p_success, size_t const num_samplings ) {
    try {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<size_t> dist(trials, p_success);
    std::map<size_t, size_t> histogram;
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    for (size_t n = 0; n < num_samplings; ++n) {
        ++histogram[dist(gen)];
    }
    for (auto p : histogram) {
        cout << std::setw(4) << p.first << ' '
                  << std::setw(10) << p.second/static_cast<double>(num_samplings) << ' '
                  << '\n';
    }
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

void binomial_rel_freq_histogram_vector(size_t const trials, double const p_success, size_t const num_samplings ) {
    try {
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::binomial_distribution<size_t> dist( trials, p_success );
    std::vector<size_t> histogram( trials+1, 0);  // create it with full size, init with 0
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    for (size_t n = 0; n < num_samplings; ++n) {
        ++histogram[ static_cast<size_t>( dist(gen) )];
    }
    size_t i = 0;
    for (auto p : histogram) {
        cout << std::setw(4) << i++ << ' '
                  << std::setw(10) << p/grostig::narrow_cast_runtime<double>(num_samplings) << ' '
                  << '\n';
    }
    //    std::copy(histogram.begin(), histogram.end(), std::ostream_iterator<int>( cout, " \n" ) );
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

class Sample_the_distribution_fo {
private: 
    Histogram_Parallel & histogram_pll;                     // this is located at begining of the class for cache performance
public:                                                     // todo: data layout for cache performance works for functions too?
    void operator()() {
        try {
        std::random_device          rd;
        std::mt19937                gen ( rd() );
        std::binomial_distribution<unsigned long>  dist { trials, p_success };  // assumed not to be thread-safe
        for (unsigned j=0; j < num_samplings; ++j) {
            ++histogram_pll[ static_cast<size_t>( dist(gen) ) + offset ];
        }
        } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
        } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
    }
    void operator()(size_t const offset_) {
        try {
        offset =                    offset_;
        std::random_device          rd;
        std::mt19937                gen ( rd() );
        std::binomial_distribution<unsigned long>  dist { trials, p_success };  // assumed not to be thread-safe
        for (unsigned j=0; j < num_samplings; ++j) {
            ++histogram_pll[ static_cast<size_t>( dist(gen) ) + offset ];
        }
        } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
        } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
    }
private:
    size_t offset {};
    size_t num_samplings {};
    size_t trials {};
    double p_success {};
    //  @OFFICER_TUBA: if only the following line is uncommented: compile error on no matching function? inititializer list?
    //                 unless it is marked static!
    //    std::random_device                  rd;     // todo: NOT USED yet!
    //    std::binomial_distribution<size_t>  dist;   // todo: NOT USED yet!
    //    std::mt19937                        gen;    // todo: NOT USED yet!
public:
    Sample_the_distribution_fo()=delete;

    explicit Sample_the_distribution_fo(Histogram_Parallel & histogram_, size_t const num_samplings_,
                               size_t const trials_, double const p_success_)
        : histogram_pll(histogram_), num_samplings(num_samplings_), trials(trials_), p_success(p_success_) {}

    explicit Sample_the_distribution_fo(Histogram_Parallel & histogram_, size_t const offset_,
                               size_t const num_samplings_, size_t const trials_, double const p_success_)
        : histogram_pll(histogram_)         /*, num_samplings(num_samplings_), offset(offset_), trials(trials_), p_success(p_success_))*/
                                            // histogram(histogram_), num_samplings(num_samplings_), offset(offset_), trials(trials_), p_success(p_success_)
    {
                                            // histogram = histogram_;      // todo: NOTE: compiler wants this initialization to be done on the "member intializer list", probably because it is a "ref", why?
        offset = offset_;
        num_samplings = num_samplings_;
        trials = trials_;
        p_success = p_success_;
        //        std::mt19937                gen ( (this->rd)() );                           // todo: NOT USED yet!
        //        std::binomial_distribution<unsigned long>  dist ( trials, p_success );      // todo: NOT USED yet!
        //                  cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    }
};

Histogram sample_the_distribution_a( size_t const num_samplings, size_t const trials, double const p_success) {
    try {
    Histogram histogram;                                    // todo: data layout?
    std::random_device          rd;
    std::mt19937                gen ( rd() );
    std::binomial_distribution<unsigned long>  dist { trials , p_success };
    for (unsigned j=0; j < num_samplings; ++j) {
        ++histogram[ static_cast<size_t>( dist(gen) ) ];
    }
    return histogram;
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
};

//void binomial_rel_freq_histogram_array_async(size_t trials, double const p_success, size_t const num_samplings ) {
//    assert( 1 <= trials && trials <= MAX_BINOMIAL_TRIALS );
//    assert( !(0.0 > p_success || p_success > 1.0) );  // between 0 and 1 ie. [0,1]
//    try {
//    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
//    //    std::random_device                          rd;                           // todo: NOT USED yet! would need to thread-safe the distribution.
//    //    std::mt19937                                gen( rd() );                  // todo: NOT USED yet!
//    //    std::binomial_distribution<size_t>          dist( trials, p_success );    // todo: NOT USED yet!
//    Histogram histogram;
//    size_t                              num_trial_values      { trials+1 };  // include the case for zero successes within the set of trials.
//    std::vector<std::future<Histogram>> futures;
//    for (size_t chunk = 0; chunk < MAX_ASYNC_FUTURES; ++chunk)
//    {
//        std::future<Histogram> future_histogram = std::async( sample_the_distribution_a, num_samplings, trials, p_success);
//        futures.push_back( future_histogram );
//    }
//    std::for_each(futures.begin(), futures.end(), // for (auto this_t : threads) {  // todo: why not?
//                  std::mem_fn( &std::thread::join ));
//    for (auto ff:futures ) {  // don't add the first chuck (ie. = 0) to itself!
//        Histogram histogram_partial = futures.pop_back().get();
//        for (size_t my_trial = 0; my_trial < num_trial_values; my_trial++) {
//            histogram[my_trial] += histogram[my_trial];
//        }
//    }
//    for (size_t i = 0; i < num_trial_values; ++i ) {
//        cout << std::setw(4) << i << ' ' << std::setw(10) << histogram[i]/static_cast<double>(num_samplings) << ' '<< '\n';
//    }
//} catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl;
//} catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
// }
//}

void binomial_rel_freq_histogram_array_thread(size_t trials, double const p_success, size_t const num_samplings ) {
    assert( 1 <= trials && trials <= MAX_BINOMIAL_TRIALS );
    assert( !(0.0 > p_success || p_success > 1.0) );  // between 0 and 1 ie. [0,1]
    try {  // sorry no indent for whole function.

    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    //    std::random_device                          rd;                           // todo: NOT USED yet! would need to thread-safe the distribution.
    //    std::mt19937                                gen( rd() );                  // todo: NOT USED yet!
    //    std::binomial_distribution<size_t>          dist( trials, p_success );    // todo: NOT USED yet!
    Histogram_Parallel              histogram;                  // this is at the top of the function for performance.
    size_t const                    chunk_sample_quantity { num_samplings / MAX_THREADS };  // todo: truncation error, off by one error
    size_t const                    num_trial_values      { trials+1 };  // include the case for zero successes within the set of trials.
    size_t const                    chunk_stride          { num_trial_values };
    size_t                          chunk_offset          { 0 };
    std::vector<std::thread>        threads;

    Sample_the_distribution_fo sample_chunk_fo2 {histogram, chunk_sample_quantity, trials, p_success};  // histogram is a reference, so we need to join within this containing function.
    for (size_t chunk = 0; chunk < MAX_THREADS; ++chunk)        //  *** Parallel Section ***
    {
        try {
            std::thread t = std::thread( sample_chunk_fo2, chunk_offset );
            threads.push_back( std::move(t) );
        } catch (...) {cerr << PGM_NAME+":fatal error: thread creation failed\n"<<endl; throw; }; // don't need thread guard due to vector. todo: even with thowing vector will be destructed and so also threads?
        chunk_offset = chunk_offset + chunk_stride;
    }                
    std::for_each(threads.begin(), threads.end(),               // *** BARRIER *** //for (auto this_t : threads) {  // todo: why not?
                  std::mem_fn( &std::thread::join ));
    chunk_offset = chunk_stride;
    for (size_t thread = 0; thread < MAX_THREADS; ++thread) {   // *** Combine Results *** // don't add the first chuck (ie. = 0) to itself!
        size_t trial_offset {0};
        for (size_t my_trial = 0; my_trial < num_trial_values; my_trial++) {
            trial_offset = chunk_offset + my_trial;
            histogram[my_trial] += histogram[trial_offset];
        }
        chunk_offset = chunk_offset + chunk_stride;
    }
    for (size_t i = 0; i < num_trial_values; ++i ) {
        cout << std::setw(4) << i << ' ' << std::setw(10) << histogram[i]/grostig::narrow_cast_runtime<double>(num_samplings) << ' '<< '\n';
    }
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

void poisson_distribution_values(double const mean, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    try {
    std::random_device rd;
    std::mt19937 gen( rd() ); // generator

    std::poisson_distribution<int> dist( mean /*4.1*/);       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto poisson = [&dist, &gen] (size_t) { return dist( gen ); };
    Eigen::RowVectorXi values = Eigen::RowVectorXi::NullaryExpr(num_values, poisson );
    std::cout << "Eigen::RowVextorXi:Poisson( mean ): (" << mean << "), "  << values << "\n";
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

void binomial_distribution_values(size_t const trials, double const p_success, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    try {
    std::random_device rd;
    std::mt19937 gen( rd() ); // generator

    std::binomial_distribution<size_t> dist( trials, p_success );       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto binomial = [&dist, &gen] (size_t) { return dist( gen ); };
    Eigen::RowVectorXi values = Eigen::RowVectorXi::NullaryExpr(num_values, binomial );
    std::cout << "Eigen::RowVextorXi:Binomial( trials, p_success ): (" << trials <<", "<< p_success <<"), "  << values << "\n";
    } catch (std::exception & e) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error exception: "<<e.what()<<endl; throw;
    } catch (...) { cerr<<PGM_NAME+":binomial_rel_freq_histogram:error unknown exception: "<<endl; throw; };
}

int main()
{
    try {
        size_t trials = 10;
        constexpr size_t  num_samplings     = 100'000'000;            // c++17
        //    binomial_rel_freq_histogram_map(    trials, 0.5, num_samplings );
        //    binomial_rel_freq_histogram_vector( trials, 0.7, num_samplings );
        binomial_rel_freq_histogram_array_thread(  trials, 1.0, num_samplings );

        poisson_distribution_values( 4.1, 100 );
        binomial_distribution_values( trials, 0.5, 100 );

        cout << "nCk(5,5): "<<boost::math::binomial_coefficient<double>(5,5)<<endl
             << "nCk(5,2): "<< boost::math::binomial_coefficient<double>(5,2)<<endl
             << "nCk(10'000,2): "<< boost::math::binomial_coefficient<double>(10000,2)<<endl;


        double mu_probability_in_bin        = 0.9;
        size_t n_samples_drawn              = 10;
        double nu_probability_of_sample     = 0.1;
        size_t k_successes = grostig::narrow_cast_runtime<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
        auto r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );

        mu_probability_in_bin               = 0.9;
        n_samples_drawn                     = 10;
        nu_probability_of_sample            = 0;
        k_successes = grostig::narrow_cast_runtime<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
        auto r2 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );
        cout << r1+r2 << endl;

        mu_probability_in_bin               = 0.9;
        n_samples_drawn                     = 10;
        nu_probability_of_sample            = 0.4;
        k_successes = grostig::narrow_cast_runtime<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
        r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );


        mu_probability_in_bin               = 0.6;
        n_samples_drawn                     = 10;
        nu_probability_of_sample            = 0;
        k_successes = grostig::narrow_cast_runtime<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
        r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );

        std::cout << "###" << std::endl;
        return 0;
    } catch (std::exception & e) {
        cerr << PGM_NAME+": exception error: "<<e.what()<<endl;
        return 1;
    } catch (...) {
        cerr << PGM_NAME+": unknown exception error.\n";
        return 2;
    }
}



