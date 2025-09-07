#include <iostream>
#include <algorithm>
#include <ranges>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>

using ld = long double;
using size_t = std::size_t;

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<ld> dist(-1., 1.);
std::uniform_real_distribution<ld> dist01(0., 1.);

struct vec2d
{
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

constexpr vec2d A_T(const ld t) {
    return vec2d{0.6 * t - 2000., 20000. - 30. * t}.normalize();
}

constexpr vec2d B(const ld t) {
    return {17188., 1750.5 - 3. * t};
}

//solve this
constexpr ld Formula1(const ld t) {
    return dot(A_T(t), B(t)) * dot(A_T(t), B(t)) - 100;
}

constexpr ld vm = 300.;
const ld theta1 = std::atan((ld) 2000. / 20000.);
const ld cos_theta1 = std::cos(theta1);
const ld sin_theta1 = std::sin(theta1);

const ld Formula2(const ld v, const ld t_explode) {
    return (20000. - 17790. - v * t_explode) / (vm * cos_theta1);
}

const ld Formula3(const ld v, const ld delta_t, const ld t_explode) {
    return (2000. - 1790. + 5 * delta_t * delta_t - 3 * t_explode) / (vm * cos_theta1 - 3.);
}

constexpr bool check_by_formula(const ld t1, const ld t2, const ld t3) {
    return t2 - t1 >= 1 && t3 - t2 >= 1;
}

namespace solve_formula1 {
    ld f(const ld t, const ld v, const ld e, const ld d) {
        ld Ax = 0.6 * t - 2000;
        ld Ay = 20000 - 30 * t;
        ld Bx = 17800 + v * e;
        ld By = 1800 - 5 * d * d - 3 * t + 3 * e;

        ld denom = std::sqrt(Ax * Ax + Ay * Ay);
        if (denom == 0) return 1e18; // 避免除零

        return (Bx * Ax + By * Ay) / denom - 10.0;
    }

    // 在区间 [L, R] 内用二分法找根（假设 f(L)*f(R) < 0）
    ld find_root(ld L, ld R, ld v, ld e, ld d) {
        for (int i = 0; i < 100; i++) { // 迭代 100 次保证精度
            ld mid = 0.5 * (L + R);
            if (f(L, v, e, d) * f(mid, v, e, d) <= 0)
                R = mid;
            else
                L = mid;
        }
        return 0.5 * (L + R);
    }

    // 主函数：返回最小实根（若不存在则返回 NaN）
    ld solve_min_root(ld v, ld e, ld d) {
        ld L = -1e6, R = 1e6; // 搜索范围可调
        ld step = 1000;       // 粗步长搜索符号变化

        std::vector<ld> roots;
        ld prev_t = L;
        ld prev_f = f(prev_t, v, e, d);

        for (ld t = L + step; t <= R; t += step) {
            ld cur_f = f(t, v, e, d);
            if (prev_f * cur_f <= 0) // 存在根
                roots.push_back(find_root(prev_t, t, v, e, d));
            prev_t = t;
            prev_f = cur_f;
        }

        std::sort(roots.begin(), roots.end());
        if (roots.empty())
            return -std::numeric_limits<ld>::infinity(); // 没有根

        const ld ans = *std::find_if(roots.begin(), roots.end(), [](ld x) { return x >= 0; });

        return ans; // 返回最小非负实根
    }
}

ld Energy(const ld v, const ld delta_t1, const ld t_explode1, const ld delta_t2, const ld t_explode2, const ld delta_t3, const ld t_explode3) {
    const ld l1 = solve_formula1::solve_min_root(v, t_explode1, delta_t1);
    const ld l2 = solve_formula1::solve_min_root(v, t_explode2, delta_t2);
    const ld l3 = solve_formula1::solve_min_root(v, t_explode3, delta_t3);
    const ld r1 = std::max(std::min(Formula2(v, t_explode1), Formula3(v, delta_t1, t_explode1)), l1);
    const ld r2 = std::max(std::min(Formula2(v, t_explode2), Formula3(v, delta_t2, t_explode2)), l2);
    const ld r3 = std::max(std::min(Formula2(v, t_explode3), Formula3(v, delta_t3, t_explode3)), l3);

    if (l1 < 0 || l2 < 0 || l3 < 0)
        return 1e18; // 无效解

    // 把区间存起来
    std::vector<std::pair<ld, ld>> intervals = {
        {l1, r1}, {l2, r2}, {l3, r3}
    };

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

void anneal() {
    ld t = 1.;
    constexpr ld cooling_rate = 0.999;
    constexpr ld min_t = 1e-12;
    ld cur_v = 105.;
    ld cur_delta_t1 = 15., cur_t_drop1 = 15.;
    ld cur_delta_t2 = 15., cur_t_drop2 = 15.;
    ld cur_delta_t3 = 15., cur_t_drop3 = 15.;

    while (t > min_t) {
        ld cur_t_explode1 = cur_t_drop1 + cur_delta_t1;
        ld cur_t_explode2 = cur_t_drop2 + cur_delta_t2;
        ld cur_t_explode3 = cur_t_drop3 + cur_delta_t3;

        ld new_v = std::clamp(cur_v + 35 * t * dist(rng), 70.L, 140.L);
        ld new_delta_t1 = std::clamp(cur_delta_t1 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop1 = std::clamp(cur_t_drop1 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_delta_t2 = std::clamp(cur_delta_t2 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop2 = std::clamp(cur_t_drop2 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_delta_t3 = std::clamp(cur_delta_t3 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop3 = std::clamp(cur_t_drop3 + 15. * t * dist(rng), 0.L, 30.L);

        ld new_t_explode1 = new_t_drop1 + new_delta_t1;
        ld new_t_explode2 = new_t_drop2 + new_delta_t2;
        ld new_t_explode3 = new_t_drop3 + new_delta_t3;

        if (!check_by_formula(new_t_drop1, new_t_drop2, new_t_drop3))
            continue;

        const ld delta = Energy(new_v, new_delta_t1, new_t_explode1, new_delta_t2, new_t_explode2, new_delta_t3, new_t_explode3) - Energy(cur_v, cur_delta_t1, cur_t_explode1, cur_delta_t2, cur_t_explode2, cur_delta_t3, cur_t_explode3);
        if (delta < 0 || std::exp(-delta / t) > dist01(rng)) {
            cur_v = new_v;
            cur_delta_t1 = new_delta_t1;
            cur_t_drop1 = new_t_drop1;
            cur_delta_t2 = new_delta_t2;
            cur_t_drop2 = new_t_drop2;
            cur_delta_t3 = new_delta_t3;
            cur_t_drop3 = new_t_drop3;
        }

        if (dist01(rng) < 0.01)
            std::cout << "Current v: " << cur_v << ", f: " << Energy(cur_v, cur_delta_t1, cur_t_explode1, cur_delta_t2, cur_t_explode2, cur_delta_t3, cur_t_explode3) << ", Temperature: " << t << "\n";

        t *= cooling_rate;
    }

    for (size_t i = 0; i < 10000; ++i) {
        ld cur_t_explode1 = cur_t_drop1 + cur_delta_t1;
        ld cur_t_explode2 = cur_t_drop2 + cur_delta_t2;
        ld cur_t_explode3 = cur_t_drop3 + cur_delta_t3;

        ld new_v = std::clamp(cur_v + 35. * t * dist(rng), 70.L, 140.L);
        ld new_delta_t1 = std::clamp(cur_delta_t1 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop1 = std::clamp(cur_t_drop1 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_delta_t2 = std::clamp(cur_delta_t2 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop2 = std::clamp(cur_t_drop2 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_delta_t3 = std::clamp(cur_delta_t3 + 15. * t * dist(rng), 0.L, 30.L);
        ld new_t_drop3 = std::clamp(cur_t_drop3 + 15. * t * dist(rng), 0.L, 30.L);

        ld new_t_explode1 = new_t_drop1 + new_delta_t1;
        ld new_t_explode2 = new_t_drop2 + new_delta_t2;
        ld new_t_explode3 = new_t_drop3 + new_delta_t3;

        if (!check_by_formula(new_t_drop1, new_t_drop2, new_t_drop3))
            continue;

        const ld delta = Energy(new_v, new_delta_t1, new_t_explode1, new_delta_t2, new_t_explode2, new_delta_t3, new_t_explode3) - Energy(cur_v, cur_delta_t1, cur_t_explode1, cur_delta_t2, cur_t_explode2, cur_delta_t3, cur_t_explode3);

        if (delta < 0) {
            cur_v = new_v;
            cur_delta_t1 = new_delta_t1;
            cur_t_drop1 = new_t_drop1;
            cur_delta_t2 = new_delta_t2;
            cur_t_drop2 = new_t_drop2;
            cur_delta_t3 = new_delta_t3;
            cur_t_drop3 = new_t_drop3;
        }

        std::cout << "Current v: " << cur_v << ", f: " << Energy(cur_v, cur_delta_t1, cur_t_explode1, cur_delta_t2, cur_t_explode2, cur_delta_t3, cur_t_explode3) << "\n";
    }

    std::cout << "cur_v: " << cur_v << ", delta_t1: " << cur_delta_t1 << ", t_drop1: " << cur_t_drop1 << ", delta_t2: " << cur_delta_t2 << ", t_drop2: " << cur_t_drop2 << ", delta_t3: " << cur_delta_t3 << ", t_drop3: " << cur_t_drop3 << "\n";
}

int main() {

    anneal();
    return 0;
}

// cur_v: 140, delta_t1: 4.25286, t_drop1: 0.360021, delta_t2: 2.67059, t_drop2: 3.30884, delta_t3: 5, t_drop3: 5