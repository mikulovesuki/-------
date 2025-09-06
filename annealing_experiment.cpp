#include <iostream>
#include <random>
#include <cmath>

using ld = long double;
using size_t = std::size_t;

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<ld> dist(-1., 1.);
std::uniform_real_distribution<ld> dist01(0., 1.);

ld f(ld x) {
    // return (x - 1) * (x - 1) * (x + 1) * (x + 1);
    // if (std::abs(x) > 100000)
        // return 1e18;
    return 2. / 3. * x * x * x * x - 1. / 6. * x * x * x - 13. / 6. * x * x + 2. / 3. * x + 2.;
}

ld anneal() {
    ld t = 100000.0;
    constexpr ld cooling_rate = 0.999;
    constexpr ld min_t = 1e-4;
    ld cur_x = 1000 + dist(rng) * 10;
    while (t > min_t) {
        ld new_x = cur_x + t * dist(rng);
        const ld delta = f(new_x) - f(cur_x);
        if (delta < 0 || std::exp(-delta / t) > dist01(rng))
            cur_x = new_x;

        if ((int) std::floor(t) % 50 == 0)
            std::cout << "Current x: " << cur_x << ", f(x): " << f(cur_x) << ", Temperature: " << t << "\n";

        t *= cooling_rate;
    }

    for (size_t i = 0; i < 10000; ++i) {
        ld new_x = cur_x + t * dist(rng);
        const ld delta = f(new_x) - f(cur_x);
        if (delta < 0) {
            cur_x = new_x;
        }

        std::cout << "Current x: " << cur_x << ", f(x): " << f(cur_x) << "\n";
    }

    return cur_x;
}

int main() {
    // anneal();
    // std::cout << f(-1) << "\n";
    // std::cout << f(0) << "\n";
    // std::cout << f(0.2) << "\n";
    // std::cout << f(1) << "\n";

    std::cout << anneal() << "\n";
    return 0;
}