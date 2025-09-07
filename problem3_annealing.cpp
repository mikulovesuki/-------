#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ranges>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <cassert>

using ld = long double;
using size_t = std::size_t;

constexpr ld ld_eq(const ld a, const ld b) {
    return std::abs(a - b) < 1e-12;
}

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<ld> dist(-1., 1.);
std::uniform_real_distribution<ld> dist01(0., 1.);

struct vec2d {
    ld x, y;

    friend constexpr ld dot(const vec2d &a, const vec2d &b);

    friend constexpr vec2d operator/(const vec2d &a, const ld b);

    constexpr ld norm() const {
        return std::hypot(x, y);
    }

    constexpr ld norm2() const {
        return dot(*this, *this);
    }

    constexpr vec2d normalize() const {
        return *this / norm();
    }
};

constexpr vec2d operator+(const vec2d &a, const vec2d &b) {
    return {a.x + b.x, a.y + b.y};
}

constexpr vec2d operator-(const vec2d &a, const vec2d &b) {
    return {a.x - b.x, a.y - b.y};
}

constexpr vec2d operator*(const ld a, const vec2d &b) {
    return {a * b.x, a * b.y};
}

constexpr vec2d operator/(const vec2d &a, const ld b) {
    return {a.x / b, a.y / b};
}

constexpr ld dot(const vec2d &a, const vec2d &b) {
    return a.x * b.x + a.y * b.y;
}

constexpr ld cross(const vec2d &a, const vec2d &b) {
    return a.x * b.y - a.y * b.x;
}

std::pair<ld, ld> quad_equation(const ld a, const ld b, const ld c) {
    const ld delta = b * b - 4 * a * c;
    if (delta < 0)
        return {NAN, NAN};
    return {(-b - std::sqrt(delta)) / (2 * a), (-b + std::sqrt(delta)) / (2 * a)};
}

const ld sqrt101 = std::sqrt(101.L);
constexpr ld g = 10.L;

constexpr ld vm = 300.L;
constexpr ld x0_m = 20000.L;
constexpr ld y0_m = 2000.L;
const ld vx_m = -300.L * 10.L / sqrt(101.L);
const ld vy_m = -300.L * 1.L / sqrt(101.L);

constexpr ld x0_fy = 17800.L;
constexpr ld y0_fy = 1800.L;

constexpr ld smoke_radius = 10.L;
constexpr ld vx_smoke = 0.L;
constexpr ld vy_smoke = -3.L;

constexpr ld x_fy(const ld v, const ld t) {
    return x0_fy + v * t;
}

const ld l_up(const ld x) {
    return x / 10 + sqrt101;
}

const ld l_down(const ld x) {
    return x / 10 - sqrt101;
}

ld calc_l(const ld v_fy, const ld t_explode, const ld t_fuze) {
    return t_explode + std::max(0.L, (y0_fy - g * t_fuze * t_fuze / 2.L - l_up(x_fy(v_fy, t_explode))) / -vy_smoke);
}

ld calc_r(const ld v_fy, const ld t_explode, const ld t_fuze, const ld l) {
    const ld ans1 = t_explode + std::max(0.L, (y0_fy - g * t_fuze * t_fuze / 2.L - l_down(x_fy(v_fy, t_explode))) / -vy_smoke);
    const ld c1 = vx_m;
    const ld c2 = x0_m - x_fy(v_fy, t_explode);
    const ld c3 = vy_m - vy_smoke;
    const ld c4 = y0_m - y0_fy + g * t_fuze * t_fuze / 2;
    const ld ans2 = quad_equation(c1 * c1 + c3 * c3, 2 * (c1 * c2 + c3 * c4), c2 * c2 + c4 * c4 - smoke_radius * smoke_radius).second;

    if (std::isnan(ans2)) {
        const ld xl_m = x0_m + l * vx_m;
        if (xl_m >= x_fy(v_fy, t_explode))
            return ans1;
        else
            return -1e18;
    }

    return ans2;
}

ld Energy(const ld v, const ld fuze_t1, const ld t_explode1, const ld fuze_t2, const ld t_explode2, const ld fuze_t3, const ld t_explode3) {
    const ld l1 = calc_l(v, t_explode1, fuze_t1);
    const ld l2 = calc_l(v, t_explode2, fuze_t2);
    const ld l3 = calc_l(v, t_explode3, fuze_t3);
    const ld r1 = calc_r(v, t_explode1, fuze_t1, l1);
    const ld r2 = calc_r(v, t_explode2, fuze_t2, l2);
    const ld r3 = calc_r(v, t_explode3, fuze_t3, l3);

    // 把区间存起来
    std::vector<std::pair<ld, ld>> intervals = {
        {l1, r1}, {l2, r2}, {l3, r3}
    };

    // 合法区间检查
    intervals.erase(std::remove_if(intervals.begin(), intervals.end(), [](const auto &interval) {
        auto [L, R] = interval;
        return L > R || std::isnan(L) || std::isnan(R) || L < 0 || R < 0;
    }), intervals.end());

    if (intervals.empty())
        return 0.L;

    // std::cout << "[" << l1 << ", " << r1 << "], [" << l2 << ", " << r2 << "], [" << l3 << ", " << r3 << "]\n";

    // 按左端点排序
    std::sort(intervals.begin(), intervals.end());

    // 合并区间
    ld total = 0.0;
    ld cur_l = intervals[0].first;
    ld cur_r = intervals[0].second;

    for (size_t i = 1; i < intervals.size(); ++i) {
        auto [L, R] = intervals[i];
        if (L > cur_r) {
            // 没有重叠，先结算前一个区间
            total += cur_r - cur_l;
            cur_l = L;
            cur_r = R;
        } else {
            // 有重叠，扩展右端点
            cur_r = std::max(cur_r, R);
        }
    }
    // 结算最后一个区间
    total += cur_r - cur_l;

    return -total; // 和原来的逻辑保持一致
}

constexpr bool check_by_formula(const ld t_drop1, const ld t_drop2, const ld t_drop3) {
    return t_drop2 - t_drop1 >= 1.L && t_drop3 - t_drop2 >= 1.L;
}

void anneal() {
    ld t = 1.L;
    constexpr ld cooling_rate = 1.L - 1e-7L;
    constexpr ld min_t = 1e-12L;
    ld cur_v = -70.L;
    ld cur_fuze_t1 = 0.L, cur_t_drop1 = 0.L;
    ld cur_fuze_t2 = 0.L, cur_t_drop2 = 0.L;
    ld cur_fuze_t3 = 0.L, cur_t_drop3 = 0.L;

    while (t > min_t) {
        ld cur_t_explode1 = cur_t_drop1 + cur_fuze_t1;
        ld cur_t_explode2 = cur_t_drop2 + cur_fuze_t2;
        ld cur_t_explode3 = cur_t_drop3 + cur_fuze_t3;

        ld new_v = std::clamp(cur_v + 5.L * t * dist(rng), -140.L, -70.L);
        ld new_fuze_t1 = std::clamp(cur_fuze_t1 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop1 = std::clamp(cur_t_drop1 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_fuze_t2 = std::clamp(cur_fuze_t2 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop2 = std::clamp(cur_t_drop2 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_fuze_t3 = std::clamp(cur_fuze_t3 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop3 = std::clamp(cur_t_drop3 + 15.L * t * dist(rng), 0.L, 30.L);

        ld new_t_explode1 = new_t_drop1 + new_fuze_t1;
        ld new_t_explode2 = new_t_drop2 + new_fuze_t2;
        ld new_t_explode3 = new_t_drop3 + new_fuze_t3;

        if (!check_by_formula(new_t_drop1, new_t_drop2, new_t_drop3))
            continue;

        const ld delta = Energy(new_v, new_fuze_t1, new_t_explode1, new_fuze_t2, new_t_explode2, new_fuze_t3, new_t_explode3) - Energy(cur_v, cur_fuze_t1, cur_t_explode1, cur_fuze_t2, cur_t_explode2, cur_fuze_t3, cur_t_explode3);
        if (delta < 0 || std::exp(-delta / t) > dist01(rng)) {
            cur_v = new_v;
            cur_fuze_t1 = new_fuze_t1;
            cur_t_drop1 = new_t_drop1;
            cur_fuze_t2 = new_fuze_t2;
            cur_t_drop2 = new_t_drop2;
            cur_fuze_t3 = new_fuze_t3;
            cur_t_drop3 = new_t_drop3;
        }

        if (dist01(rng) < 1e-4L) {
            std::cout << std::fixed << std::setprecision(6) << "Current v: " << cur_v << ", f: " << Energy(cur_v, cur_fuze_t1, cur_t_explode1, cur_fuze_t2, cur_t_explode2, cur_fuze_t3, cur_t_explode3);
            std::cout << std::defaultfloat << ", Temperature: " << t << "\n";
        }

        t *= cooling_rate;
    }

    for (size_t i = 0; i < 10000; ++i) {
        ld cur_t_explode1 = cur_t_drop1 + cur_fuze_t1;
        ld cur_t_explode2 = cur_t_drop2 + cur_fuze_t2;
        ld cur_t_explode3 = cur_t_drop3 + cur_fuze_t3;

        ld new_v = std::clamp(cur_v + 35.L * t * dist(rng), -140.L, -70.L);
        ld new_fuze_t1 = std::clamp(cur_fuze_t1 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop1 = std::clamp(cur_t_drop1 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_fuze_t2 = std::clamp(cur_fuze_t2 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop2 = std::clamp(cur_t_drop2 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_fuze_t3 = std::clamp(cur_fuze_t3 + 15.L * t * dist(rng), 0.L, 30.L);
        ld new_t_drop3 = std::clamp(cur_t_drop3 + 15.L * t * dist(rng), 0.L, 30.L);

        ld new_t_explode1 = new_t_drop1 + new_fuze_t1;
        ld new_t_explode2 = new_t_drop2 + new_fuze_t2;
        ld new_t_explode3 = new_t_drop3 + new_fuze_t3;

        if (!check_by_formula(new_t_drop1, new_t_drop2, new_t_drop3))
            continue;

        const ld delta = Energy(new_v, new_fuze_t1, new_t_explode1, new_fuze_t2, new_t_explode2, new_fuze_t3, new_t_explode3) - Energy(cur_v, cur_fuze_t1, cur_t_explode1, cur_fuze_t2, cur_t_explode2, cur_fuze_t3, cur_t_explode3);

        if (delta < 0) {
            cur_v = new_v;
            cur_fuze_t1 = new_fuze_t1;
            cur_t_drop1 = new_t_drop1;
            cur_fuze_t2 = new_fuze_t2;
            cur_t_drop2 = new_t_drop2;
            cur_fuze_t3 = new_fuze_t3;
            cur_t_drop3 = new_t_drop3;
        }

        std::cout << std::fixed << std::setprecision(6) << "Current v: " << cur_v << ", f: " << Energy(cur_v, cur_fuze_t1, cur_t_explode1, cur_fuze_t2, cur_t_explode2, cur_fuze_t3, cur_t_explode3) << "\n";
    }

    std::cout << std::fixed << std::setprecision(6) << "cur_v: " << cur_v << ", fuze_t1: " << cur_fuze_t1 << ", t_drop1: " << cur_t_drop1 << ", fuze_t2: " << cur_fuze_t2 << ", t_drop2: " << cur_t_drop2 << ", fuze_t3: " << cur_fuze_t3 << ", t_drop3: " << cur_t_drop3 << "\n" << "f: " << Energy(cur_v, cur_fuze_t1, cur_t_drop1 + cur_fuze_t1, cur_fuze_t2, cur_t_drop2 + cur_fuze_t2, cur_fuze_t3, cur_t_drop3 + cur_fuze_t3) << "\n";
    std::cout << "t1: " << calc_r(cur_v, cur_t_drop1 + cur_fuze_t1, cur_fuze_t1, calc_l(cur_v, cur_t_drop1 + cur_fuze_t1, cur_fuze_t1)) - calc_l(cur_v, cur_t_drop1 + cur_fuze_t1, cur_fuze_t1) << "\n";
    std::cout << "t2: " << calc_r(cur_v, cur_t_drop2 + cur_fuze_t2, cur_fuze_t2, calc_l(cur_v, cur_t_drop2 + cur_fuze_t2, cur_fuze_t2)) - calc_l(cur_v, cur_t_drop2 + cur_fuze_t2, cur_fuze_t2) << "\n";
    std::cout << "t3: " << calc_r(cur_v, cur_t_drop3 + cur_fuze_t3, cur_fuze_t3, calc_l(cur_v, cur_t_drop3 + cur_fuze_t3, cur_fuze_t3)) - calc_l(cur_v, cur_t_drop3 + cur_fuze_t3, cur_fuze_t3) << "\n";
}

int main() {
    anneal();
    return 0;
}

// cur_v: -99.078733, fuze_t1: 2.714645, t_drop1: 0.000000, fuze_t2: 3.587747, t_drop2: 2.311771, fuze_t3: 4.883583, t_drop3: 6.147723
// f: -15.016578
// t1: 6.699917
// t2: 6.699917
// t3: 6.699917