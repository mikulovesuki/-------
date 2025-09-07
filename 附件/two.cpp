// sa_three_smokes_fixed.cpp
#include <bits/stdc++.h>
using namespace std;

/* ------------------ 常量（题目给定） ------------------ */
const double g = 9.8;
const double R_smoke = 10.0;
const double v_smoke = 3.0;
const double smoke_active = 20.0;
const double vM = 300.0;
const array<double,3> M0 = {20000.0, 0.0, 2000.0};
const array<double,3> fake_target = {0.0, 0.0, 0.0};
const array<double,3> U0 = {17800.0, 0.0, 1800.0};
const array<double,3> C_true = {0.0, 200.0, 0.0};

/* ------------------ 向量工具 ------------------ */
using Vec3 = array<double,3>;
inline Vec3 addv(const Vec3 &a,const Vec3 &b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
inline Vec3 subv(const Vec3 &a,const Vec3 &b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline Vec3 scalv(const Vec3 &a,double s){return {a[0]*s,a[1]*s,a[2]*s};}
inline double dotp(const Vec3 &a,const Vec3 &b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline double norm2(const Vec3 &a){return dotp(a,a);}
inline double norm(const Vec3 &a){return sqrt(norm2(a));}

/* 计算导弹速度向量（恒向量） */
Vec3 missile_velocity_vector() {
    Vec3 d = {fake_target[0]-M0[0], fake_target[1]-M0[1], fake_target[2]-M0[2]};
    double n = norm(d);
    return scalv(d, vM / n);
}

/* 导弹位置 */
Vec3 missile_pos(double t, const Vec3 &dirM) {
    return addv(M0, scalv(dirM, t));
}

/* 无人机位置（等高直线） */
Vec3 UAV_pos(double theta, double v, double t) {
    return { U0[0] + v * cos(theta) * t,
             U0[1] + v * sin(theta) * t,
             U0[2] };
}

/* 起爆点 E（基于投放时刻 t_drop 与起爆延迟 t_expl） */
Vec3 compute_E_point(double theta, double v, double t_drop, double t_expl) {
    Vec3 B0 = UAV_pos(theta, v, t_drop);
    double Ex = B0[0];
    double Ey = B0[1];
    double Ez = B0[2] - 0.5 * g * t_expl * t_expl;
    return {Ex, Ey, Ez};
}

/* 云团当前中心（考虑下沉） */
Vec3 smoke_center_at(const Vec3 &E, double t_explode, double t_now) {
    double dt = t_now - t_explode;
    if (dt < 0.0) dt = 0.0;
    return { E[0], E[1], E[2] - v_smoke * dt };
}

/* 判断时刻 t_now 是否被某枚烟幕遮蔽（线段 M--C_true 与球相交，且烟幕有效） */
bool isOccludedAtTime(const Vec3 &Mpos, const Vec3 &Gpos, const Vec3 &E, double t_explode, double t_now) {
    const double EPS = 1e-9;
    if (t_now < t_explode - EPS) return false;
    if (t_now > t_explode + smoke_active + EPS) return false;
    Vec3 C = smoke_center_at(E, t_explode, t_now);
    Vec3 d = { Gpos[0] - Mpos[0], Gpos[1] - Mpos[1], Gpos[2] - Mpos[2] };
    double A = norm2(d);
    if (A < EPS) {
        double dist2 = norm2({Mpos[0]-C[0], Mpos[1]-C[1], Mpos[2]-C[2]});
        return dist2 <= R_smoke*R_smoke + EPS;
    }
    Vec3 MC = { Mpos[0]-C[0], Mpos[1]-C[1], Mpos[2]-C[2] };
    double B = 2.0 * (d[0]*MC[0] + d[1]*MC[1] + d[2]*MC[2]);
    double C0 = norm2(MC) - R_smoke*R_smoke;
    double discr = B*B - 4.0*A*C0;
    if (discr < -EPS) return false;
    if (discr < 0) discr = 0;
    double sqrtD = sqrt(discr);
    double u1 = (-B - sqrtD) / (2.0*A);
    double u2 = (-B + sqrtD) / (2.0*A);
    const double LO = -1e-9, HI = 1.0 + 1e-9;
    return ( (u1>=LO && u1<=HI) || (u2>=LO && u2<=HI) );
}

/* 计算三枚烟幕在时间轴上的并集遮蔽时长（数值积分） */
double union_occlusion_time(double theta, double v,
                            const array<double,3> &drops, const array<double,3> &expls,
                            double tmax = 100.0, double dt = 0.1)
{
    // 计算 E 点与 T_explode
    array<Vec3,3> Es;
    array<double,3> T_ex;
    for (int i=0;i<3;i++){
        Es[i] = compute_E_point(theta, v, drops[i], expls[i]);
        T_ex[i] = drops[i] + expls[i];
    }
    Vec3 dirM = missile_velocity_vector();
    int nsteps = (int)(tmax/dt) + 1;
    vector<char> oc(nsteps, 0);
    for (int k=0;k<nsteps;k++){
        double t = k*dt;
        Vec3 Mpos = missile_pos(t, dirM);
        for (int i=0;i<3;i++){
            if (isOccludedAtTime(Mpos, C_true, Es[i], T_ex[i], t)) { oc[k]=1; break; }
        }
    }
    int cnt = 0; for (char b: oc) if (b) cnt++;
    return cnt * dt;
}

/* 目标函数（带罚项） */
double objective_with_penalty(double theta, double v,
                              const array<double,3> &drops, const array<double,3> &expls)
{
    double penalty = 0.0;
    // enforce spacing between drops >=1s (on sorted drops)
    vector<int> sorted_idx = {0,1,2};
    sort(sorted_idx.begin(), sorted_idx.end(),
         [&](int a,int b){ return drops[a] < drops[b]; });
    for (int k=0;k<2;k++){
        double gap = drops[sorted_idx[k+1]] - drops[sorted_idx[k]];
        if (gap < 1.0) penalty += 100.0 * (1.0 - gap);
    }
    if (v < 70.0 || v > 140.0) penalty += 50.0 * (fabs(min(0.0, v-70.0)) + fabs(max(0.0, v-140.0)));
    for (int i=0;i<3;i++){
        if (drops[i] < 0.0) penalty += 100.0 * (-drops[i]);
        if (expls[i] < 0.0) penalty += 100.0 * (-expls[i]);
    }
    double oc = union_occlusion_time(theta, v, drops, expls, 100.0, 0.1);
    return oc - penalty;
}

/* 随机生成初始解与邻域函数 */
double randd(double a,double b){ return a + (b-a) * (rand() / (double)RAND_MAX); }

struct State {
    double theta;
    double v;
    array<double,3> drops;
    array<double,3> expls;
};

State random_initial_state() {
    State s;
    s.theta = randd(-M_PI, M_PI);
    s.v = randd(70.0,140.0);
    vector<double> d(3);
    for (int i=0;i<3;i++) d[i] = randd(0.0,40.0);
    sort(d.begin(), d.end());
    for (int i=0;i<3;i++) s.drops[i] = d[i];
    for (int i=0;i<3;i++) s.expls[i] = randd(0.1,8.0);
    return s;
}

State neighbor_state(const State &s) {
    State ns = s;
    ns.theta += randd(-0.3, 0.3);
    ns.v += randd(-8.0, 8.0);
    for (int i=0;i<3;i++){
        ns.drops[i] += randd(-1.0, 1.0);
        ns.expls[i] += randd(-0.5, 0.5);
        if (ns.drops[i] < 0.0) ns.drops[i] = 0.0;
        if (ns.expls[i] < 0.0) ns.expls[i] = 0.0;
    }
    if (ns.v < 50.0) ns.v = 50.0;
    if (ns.v > 160.0) ns.v = 160.0;
    return ns;
}

/* 模拟退火流程 */
State simulated_annealing(int iterations=4000, double T0=1.0, double alpha=0.998) {
    State cur = random_initial_state();
    double cur_val = objective_with_penalty(cur.theta, cur.v, cur.drops, cur.expls);
    State best = cur; double best_val = cur_val;
    double T = T0;
    for (int it=0; it<iterations; ++it) {
        State ns = neighbor_state(cur);
        double nv = objective_with_penalty(ns.theta, ns.v, ns.drops, ns.expls);
        double delta = nv - cur_val;
        if (delta > 0 || exp(delta / (T + 1e-12)) > (rand() / (double)RAND_MAX)) {
            cur = ns; cur_val = nv;
            if (cur_val > best_val) { best = cur; best_val = cur_val; }
        }
        T *= alpha;
    }
    // print best objective found
    cout << fixed << setprecision(6);
    cout << "Best objective (union occlusion seconds - penalty) = " << best_val << endl;
    double true_oc = union_occlusion_time(best.theta, best.v, best.drops, best.expls, 100.0, 0.1);
    cout << "Recomputed true union occlusion (s): " << true_oc << endl;
    return best;
}

/* ------------------ main ------------------ */
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    srand((unsigned)time(nullptr));

    State ans = simulated_annealing(4000, 1.0, 0.998);

    cout << fixed << setprecision(6);
    cout << "Final strategy (theta [rad], theta[deg], v):\n";
    cout << ans.theta << " , " << ans.theta*180.0/M_PI << " , " << ans.v << endl;
    for (int i=0;i<3;i++){
        double Tex = ans.drops[i] + ans.expls[i];
        Vec3 E = compute_E_point(ans.theta, ans.v, ans.drops[i], ans.expls[i]);
        cout << "smoke " << (i+1) << " : drop=" << ans.drops[i]
             << " expl_delay=" << ans.expls[i]
             << " explode_time=" << Tex
             << " E=("<<E[0]<<","<<E[1]<<","<<E[2]<<")\n";
    }
    // recompute and print individual contributions
    double dt = 0.1, tmax = 100.0;
    int nsteps = (int)(tmax/dt)+1;
    vector<char> unionMask(nsteps,0);
    vector<vector<char>> per(3, vector<char>(nsteps,0));
    array<Vec3,3> Es;
    array<double,3> T_ex;
    for (int i=0;i<3;i++){ Es[i]=compute_E_point(ans.theta, ans.v, ans.drops[i], ans.expls[i]); T_ex[i]=ans.drops[i]+ans.expls[i]; }
    Vec3 dirM = missile_velocity_vector();
    for (int k=0;k<nsteps;k++){
        double t = k*dt;
        Vec3 Mpos = missile_pos(t, dirM);
        for (int i=0;i<3;i++){
            if (isOccludedAtTime(Mpos, C_true, Es[i], T_ex[i], t)){ per[i][k]=1; unionMask[k]=1; }
        }
    }
    double unionSec = 0.0;
    vector<double> indiv(3,0.0);
    for (int k=0;k<nsteps;k++){
        if (unionMask[k]) unionSec += dt;
        for (int i=0;i<3;i++) if (per[i][k]) indiv[i] += dt;
    }
    cout << "Total union occlusion (s) = " << unionSec << endl;
    for (int i=0;i<3;i++) cout << "smoke " << (i+1) << " individual occlusion(s) = " << indiv[i] << endl;

    return 0;
}
