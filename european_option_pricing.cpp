#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>

// Function to generate random numbers from a normal distribution
double generate_normal_random() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<> d(0,1);
    return d(gen);
}

// Geometric Brownian Motion with varying interest rate
std::vector<double> GBM_varying_interest(int days, double initial_price, double r_0, double alpha, double sigmar, double std) {
    std::vector<double> St(days);
    std::vector<double> r(days);
    St[0] = initial_price;
    r[0] = r_0;

    for (int i = 1; i < days; ++i) {
        r[i] = r[i-1] * std::exp(-alpha) + r_0 * (1 - std::exp(-alpha)) +
               sigmar * std::sqrt((1 - std::exp(-2 * alpha)) / (2 * alpha)) * generate_normal_random();
        
        St[i] = St[i-1] * std::exp((r[i] - 0.5 * std::pow(std, 2)) + std * generate_normal_random());
    }

    return St;
}

// Geometric Brownian Motion with varying volatility
std::vector<double> GBM_varying_volatility(int days, double initial_price, double r, double sigma_0, double a_sigma, double b_sigma) {
    std::vector<double> St(days);
    std::vector<double> sigma(days);
    St[0] = initial_price;
    sigma[0] = sigma_0;

    for (int j = 1; j < days; ++j) {
        sigma[j] = sigma_0 + a_sigma * ((St[j-1] / initial_price) - 1) +
                   b_sigma * std::pow((St[j-1] / initial_price) - 1, 2);
        St[j] = St[j-1] * std::exp((r - 0.5 * std::pow(sigma[j], 2))  + 
                                 sigma[j]  * generate_normal_random());
    }
    return St;
}

// Monte Carlo simulation for European options
double monte_carlo_european_option(int num_simulations, std::string option_type, double S0, double K, int days, double r, std::string GeometricBrownianMotion_type, double initial_price, double r_0, double alpha, double sigmar, double std, double sigma_0, double a_sigma, double b_sigma) {
    std::vector<double> payoffs(num_simulations);

    for (int i = 0; i < num_simulations; ++i) {
        std::vector<double> St_path;
        if (GeometricBrownianMotion_type == "varying_interest") {
            St_path = GBM_varying_interest(days, initial_price, r_0, alpha, sigmar, std);
        } else if (GeometricBrownianMotion_type == "varying_volatility") {
            St_path = GBM_varying_volatility(days, initial_price, r, sigma_0, a_sigma, b_sigma);
        }
        double St_final = St_path.back();
        
        if (option_type == "call") {
            payoffs[i] = std::max(St_final - K, 0.0);
        } else if (option_type == "put") {
            payoffs[i] = std::max(K - St_final, 0.0);
        }
    }

    double average_discounted_payoff = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / num_simulations;
    return average_discounted_payoff * std::exp(-r * (days / 365.0));
}

int main() {
    // Parameters for simulations
    double initial_price = 100;
    double r_0 = 0.05;
    double alpha = 0.2;
    double sigmar = 0.03;
    double std = 0.2;
    double sigma_0 = 0.2;
    double a_sigma = 0.1;
    double b_sigma = 0.05;

    // Running test simulations
    double payoff_constant_volatility = monte_carlo_european_option(10000000, "call", 100, 105, 90, 0.05, "varying_interest", initial_price, r_0, alpha, sigmar, std, sigma_0, a_sigma, b_sigma);
    double payoff_varying_volatility = monte_carlo_european_option(10000000, "call", 100, 105, 90, 0.05, "varying_volatility", initial_price, r_0, alpha, sigmar, std, sigma_0, a_sigma, b_sigma);

    std::cout << "Payoff with constant volatility: " << payoff_constant_volatility << "\n";
    std::cout << "Payoff with varying volatility: " << payoff_varying_volatility << std::endl;

    return 0;
}