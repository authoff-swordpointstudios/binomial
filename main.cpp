#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <unordered_map>
#include <iterator>
#include <random>

#include <boost/math/special_functions/binomial.hpp>
#include <tr1/cmath>  // to pick up beta()
#include <eigen3/Eigen/Dense>

using std::cout;
using std::cin;
using std::endl;

//double binomial_coefficient_nCk(int const n, int const k) {  // THIS IS partly WRONG, use boost.
//    if (n == 0 || n == k) return 1;
//    return 1/ ((n+1) * std::tr1::beta( n-k+1, k+1 ));
//}

double binomial_probability_mass_function( size_t const k_successes, size_t const n_trials, double const probability_of_success) {
    auto bc =           boost::math::binomial_coefficient<double>( 5, 2);
    auto p_k =          std::pow(probability_of_success, k_successes);
    auto p_complement = std::pow(1 - probability_of_success, n_trials - k_successes);
    return bc * p_k * p_complement;
}

void histogram_binomial_map(size_t const trials, double const p_success ) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<size_t> dist(trials, p_success);
    std::map<size_t, size_t> histogram;
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;

    for (int n = 0; n < 1'000'000; ++n) {
        ++histogram[dist(gen)];
    }
    for (auto p : histogram) {
        cout << std::setw(4) << p.first << ' '
                  << std::setw(10) << p.second << ' '
                  << '\n';
    }
}

void histogram_binomial_vector(size_t const trials, double const p_success ) {
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::binomial_distribution<size_t> dist( trials, p_success );
    std::vector<int> histogram( static_cast<size_t>( trials + 1 ) , 0);
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;

    double probes = 1'000'000;
    for (long n = 0; n < probes; ++n) {
        ++histogram[ static_cast<size_t>( dist(gen) )];
    }
    int i = 0;
    for (auto p : histogram) {
        cout << std::setw(4) << i++ << ' '
                  << std::setw(10) << p/probes << ' '
                  << '\n';
    }
//    std::copy(histogram.begin(), histogram.end(), std::ostream_iterator<int>( cout, " \n" ) );
}

void histogram_binomial_array(size_t trials, double const p_success ) {
    constexpr int MAX_TRIALS = 1000;
    assert(1 <= trials && trials <= MAX_TRIALS);
    assert( !(0.0 > p_success || p_success > 1.0) );
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::binomial_distribution<size_t> dist( trials, p_success );
    std::array<int, MAX_TRIALS> histogram;
    cout << "Probability Mass/Density Function Graph - Binomial Distribution with trials, likelyhood: " << trials << ", " << p_success << std::endl;

    double probes = 1'000'000;
    for (long n = 0; n < probes; ++n) {
        ++histogram[ static_cast<size_t>( dist(gen) )];
    }
    for (size_t i = 0; i < trials+1; ++i ) {
        cout << std::setw(4) << i << ' '
                  << std::setw(10) << histogram[i]/probes << ' '
                  << '\n';
    }
}

void distribution_values_binomial(double const mean, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    std::random_device rd;
    std::mt19937 gen( rd() ); // generator

    std::poisson_distribution<int> dist( mean /*4.1*/);       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto poisson = [&dist, &gen] (int) { return dist( gen ); };
    Eigen::RowVectorXi v = Eigen::RowVectorXi::NullaryExpr(num_values, poisson );
    std::cout << "Eigen::RowVextorXi:Poisson( mean ): (" << mean << "), "  << v << "\n";
}

void distribution_values_binomial(size_t const trials, double const p_success, long int const num_values = 100 ) { // creates a vector loaded with values from the distribution during its initialization.
    std::random_device rd;
    std::mt19937 gen(rd()); // generator

    std::binomial_distribution<size_t> dist( trials, p_success );       // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto binomial = [&dist, &gen] (size_t) { return dist( gen ); };
    Eigen::RowVectorXi v = Eigen::RowVectorXi::NullaryExpr(num_values, binomial );
    std::cout << "Eigen::RowVextorXi:Binomial( trials, p_success ): (" << trials <<", "<< p_success <<"), "  << v << "\n";
}

double binomial_theorem(double const mu_probability_in_bin, size_t const n_samples_drawn, size_t const k_successes)
{
    auto r1 = binomial_probability_mass_function(k_successes, n_samples_drawn, mu_probability_in_bin);
    cout << "binomial_probability_mass_function(int k_successes, int n_trials, double p): " << k_successes << ", " << n_samples_drawn << ", " << mu_probability_in_bin << endl;
    cout << r1 << endl;
    return r1;
}

int main()
{
    cout << "nCk(5,5): "<<boost::math::binomial_coefficient<double>(5,5)<<endl
         << "nCk(5,2): "<< boost::math::binomial_coefficient<double>(5,2)<<endl
         << "nCk(10'000,2): "<< boost::math::binomial_coefficient<double>(10'000,2)<<endl;
    size_t trials = 10;
    histogram_binomial_array( trials, 0.5 );
    histogram_binomial_vector( trials, 0.7 );
    histogram_binomial_vector( trials, 1.0 );
    distribution_values_binomial( 4.1, 100 );
    distribution_values_binomial( trials, 0.5, 100 );

    double mu_probability_in_bin        = 0.9;
    size_t n_samples_drawn                 = 10;
    double nu_probability_of_sample        = 0.1;
    size_t k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    auto r1 = binomial_theorem( mu_probability_in_bin, n_samples_drawn, k_successes);


    mu_probability_in_bin               = 0.9;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));
    auto r2 = binomial_theorem( mu_probability_in_bin, n_samples_drawn, k_successes);
    cout << r1+r2 << endl;

    mu_probability_in_bin               = 0.9;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0.4;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));

    mu_probability_in_bin               = 0.9;
    r1 = binomial_theorem( mu_probability_in_bin, n_samples_drawn, k_successes);


    mu_probability_in_bin               = 0.6;
    n_samples_drawn                     = 10;
    nu_probability_of_sample            = 0;
    k_successes = static_cast<size_t>( lround( nu_probability_of_sample * n_samples_drawn ));

    r2 = binomial_theorem( mu_probability_in_bin, n_samples_drawn, k_successes);


    std::cout << "###" << std::endl;
    return 0;
}

/* std::cout << "Pascal's triangle:\n";
for(int n = 1; n < 10; ++n) {
    std::cout << std::string(20-n*2, ' ');
    for(int k = 1; k < n; ++k)
        std::cout << std::setw(3) << binomial_coefficient_nCk(n,k) << ' ';
    std::cout << '\n';
} */


