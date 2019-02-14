#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <unordered_map>
#include <iterator>
#include <random>
#include <thread>
#include <algorithm>
#include <functional>

#include <boost/math/special_functions/binomial.hpp>
//#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <tr1/cmath>  // to pick up beta()
#include <eigen3/Eigen/Dense>
using std::cout;
using std::cin;
using std::endl;
constexpr int    NUM_CPU_CORES   = 8;
constexpr int    MAX_THREADS     = NUM_CPU_CORES;
constexpr size_t MAX_BINOMIAL_TRIALS = 1000;
constexpr size_t TRIALS_ARRAY_SZ = (MAX_BINOMIAL_TRIALS+1)*MAX_THREADS;

void boost_log_init() {  // BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
        boost::log::core::get()->set_filter( boost::log::trivial::severity >= boost::log::trivial::info);
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

double binomial_distribution_probability_mass_function( size_t const n_trials, size_t const k_successes, double const probability_of_success) {
    auto bc =           boost::math::binomial_coefficient<double>( 5, 2);
    auto p_k =          std::pow(probability_of_success, k_successes);
    auto p_complement = std::pow(1 - probability_of_success, n_trials - k_successes);
    return bc * p_k * p_complement;
}

double print_binomial_dist_PMF(size_t const n_samples_drawn, size_t const k_successes, double const mu_probability_in_bin)
{
    auto r1 = binomial_distribution_probability_mass_function(n_samples_drawn, k_successes,  mu_probability_in_bin);
    cout << "binomial_probability_mass_function(n_trials,k_successes,p): " << n_samples_drawn << ", " << k_successes << ", " << mu_probability_in_bin << endl;
    cout << r1 << endl;
    return r1;
}

void binomial_rel_freq_histogram_map(size_t const trials, double const p_success, size_t const num_samplings ) {
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
}

void binomial_rel_freq_histogram_vector(size_t const trials, double const p_success, size_t const num_samplings ) {
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::binomial_distribution<size_t> dist( trials, p_success );
    std::vector<size_t> histogram( static_cast<size_t>( trials + 1 ) , 0);  // create it with full size, init with 0
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    for (size_t n = 0; n < num_samplings; ++n) {
        ++histogram[ static_cast<size_t>( dist(gen) )];
    }
    size_t i = 0;
    for (auto p : histogram) {
        cout << std::setw(4) << i++ << ' '
                  << std::setw(10) << p/static_cast<double>(num_samplings) << ' '
                  << '\n';
    }
    //    std::copy(histogram.begin(), histogram.end(), std::ostream_iterator<int>( cout, " \n" ) );
}

class Sample_the_distribution_fo {
private: 
    using Hist = std::array< size_t, TRIALS_ARRAY_SZ >;     // todo : gets quite large...
    Hist & histogram;                                       // this is located at begining of the class for cache performance
public:                                                     // todo: data layout for cache performance works for functions too?
    void operator()() {
        std::random_device          rd;
        std::mt19937                gen ( rd() );
        std::binomial_distribution<unsigned long>  dist { trials , p_success };
        for (unsigned j=0; j < num_samplings; ++j) {
            ++histogram[ static_cast<size_t>( dist(gen) ) + offset ];
        }
    }
private:
    size_t offset;
    size_t num_samplings;
    size_t trials;
    double p_success;
    //  @OFFICER_TUBA: if only the following line is uncommented: compile error on no matching function? inititializer list?
    //                 unless it is marked static!
    //    std::random_device                  rd;     // todo: NOT USED yet!
    //    std::binomial_distribution<size_t>  dist;   // todo: NOT USED yet!
    //    std::mt19937                        gen;    // todo: NOT USED yet!
public:
    Sample_the_distribution_fo()=delete;                // todo: is this a good idea?
    Sample_the_distribution_fo(Hist & histogram_, size_t offset_, size_t const num_samplings_, size_t const trials_, double const p_success_)
        : histogram(histogram_)/*, num_samplings(num_samplings_), offset(offset_), trials(trials_), p_success(p_success_))*/
        // histogram(histogram_), num_samplings(num_samplings_), offset(offset_), trials(trials_), p_success(p_success_)
    {
        // histogram = histogram_;      // todo: NOTE: compiler wants this initialization to be done on the "member intializer list", probably because it is a "ref"
        offset = offset_;
        num_samplings = num_samplings_;
        trials = trials_;
        p_success = p_success_;
        //        std::mt19937                gen ( (this->rd)() );                           // todo: NOT USED yet!
        //        std::binomial_distribution<unsigned long>  dist ( trials, p_success );      // todo: NOT USED yet!
        //                  cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    }
};

void binomial_rel_freq_histogram_array(size_t trials, double const p_success, size_t const num_samplings ) {
    assert( 1 <= trials && trials <= MAX_BINOMIAL_TRIALS );
    assert( !(0.0 > p_success || p_success > 1.0) );  // between 0 and 1 ie. [0,1]
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;
    //    std::random_device                          rd;                           // todo: NOT USED yet! would need to thread-safe the distribution.
    //    std::mt19937                                gen( rd() );                  // todo: NOT USED yet!
    //    std::binomial_distribution<size_t>          dist( trials, p_success );    // todo: NOT USED yet!
    std::array<size_t, TRIALS_ARRAY_SZ> histogram;
    size_t                              chunk_sample_quantity { num_samplings / MAX_THREADS };  // todo: truncation error, off by one error
    size_t                              num_trial_values      { trials+1 };  // include the case for zero successes within the set of trials.
    size_t                              chunk_offset          { 0 };
    size_t                              chunk_stride          { num_trial_values };
    std::vector<std::thread>            threads;

    for (size_t chunk = 0; chunk < MAX_THREADS; ++chunk)
    {
        Sample_the_distribution_fo sample_chunk_fo2{ histogram, chunk_offset, chunk_sample_quantity, trials, p_success };
                  //        std::thread t = std::thread( sample_chunk_fo2 );
                  //        std::thread t = std::thread( sample_chunk_fo, histogram, offset, num_samplings/(chunk+1), dist, gen );
                  //        threads.push_back( t );
        threads.push_back( std::thread( (sample_chunk_fo2) ) );  // todo:  why can't I put t in here?  get a deleted function template error
                  //        thread_guard guarded_thread(t);  // todo: yes it is guarded.  how do I get return value back from thread? you don't with threads.
        chunk_offset = chunk_offset + chunk_stride;
    }                
    std::for_each(threads.begin(), threads.end(), // for (auto this_t : threads) {  // todo: why not?
                  std::mem_fn( &std::thread::join ));
    for (size_t i = 0; i < /*trials+1*/ 110; ++i ) {
        cout << std::setw(4) << i << ' '
                  << std::setw(10) << histogram[i]/*/num_samplings*/ << ' '
                  << '\n';
    }
    chunk_offset = chunk_stride;
    for (size_t thread = 0; thread < MAX_THREADS; ++thread) {  // don't add the first chuck (ie. = 0) to itself!
        size_t trial_offset {0};
        for (size_t my_trial = 0; my_trial < num_trial_values; my_trial++) {
            trial_offset = chunk_offset + my_trial;
            histogram[my_trial] += histogram[trial_offset];
        }
        chunk_offset = chunk_offset + chunk_stride;
    }
    for (size_t i = 0; i < /*trials+1*/ 110; ++i ) {
        cout << std::setw(4) << i << ' '
                  << std::setw(10) << histogram[i]/static_cast<double>(num_samplings) << ' '
                  << '\n';
    }
}

void poisson_distribution_values(double const mean, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    std::random_device rd;
    std::mt19937 gen( rd() ); // generator

    std::poisson_distribution<int> dist( mean /*4.1*/);       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto poisson = [&dist, &gen] (size_t) { return dist( gen ); };
    Eigen::RowVectorXi values = Eigen::RowVectorXi::NullaryExpr(num_values, poisson );
    std::cout << "Eigen::RowVextorXi:Poisson( mean ): (" << mean << "), "  << values << "\n";
}

void binomial_distribution_values(size_t const trials, double const p_success, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    std::random_device rd;
    std::mt19937 gen( rd() ); // generator

    std::binomial_distribution<size_t> dist( trials, p_success );       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto binomial = [&dist, &gen] (size_t) { return dist( gen ); };
    Eigen::RowVectorXi values = Eigen::RowVectorXi::NullaryExpr(num_values, binomial );
    std::cout << "Eigen::RowVextorXi:Binomial( trials, p_success ): (" << trials <<", "<< p_success <<"), "  << values << "\n";
}

int main()
{
    size_t trials = 10;
    constexpr size_t  num_samplings     = 100'000'000;            // c++17
//    binomial_rel_freq_histogram_map(    trials, 0.5, num_samplings );
//    binomial_rel_freq_histogram_vector( trials, 0.7, num_samplings );
    binomial_rel_freq_histogram_array(  trials, 1.0, num_samplings );

    poisson_distribution_values( 4.1, 100 );
    binomial_distribution_values( trials, 0.5, 100 );

    cout << "nCk(5,5): "<<boost::math::binomial_coefficient<double>(5,5)<<endl
         << "nCk(5,2): "<< boost::math::binomial_coefficient<double>(5,2)<<endl
         << "nCk(10'000,2): "<< boost::math::binomial_coefficient<double>(10000,2)<<endl;


    double mu_probability_in_bin        = 0.9;
    size_t n_samples_drawn              = 10;
    double nu_probability_of_sample     = 0.1;
    size_t k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    auto r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );


    mu_probability_in_bin               = 0.9;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    auto r2 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );
    cout << r1+r2 << endl;

    mu_probability_in_bin               = 0.9;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0.4;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );


    mu_probability_in_bin               = 0.6;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    r1 = print_binomial_dist_PMF( n_samples_drawn, k_successes, mu_probability_in_bin );

    std::cout << "###" << std::endl;
    return 0;
}



