#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using std::complex;
using std::cout;
using std::endl;
using std::string;
using std::vector;

constexpr double PI = 3.14159265358979323846;

struct Branch {
    int i;
    int j;
    double r;
    double x;
};

struct NRSolution {
    bool converged = false;
    int iterations = 0;
    vector<double> Vm;
    vector<double> Va;
    vector<double> Pcalc;
    vector<double> Qcalc;
    vector<int> activePv;
    vector<int> activePq;
    vector<int> qLimited;
};

struct Frame {
    double hour = 0.0;
    vector<double> Vm;
    vector<double> VaDeg;
    vector<double> Pload;
    vector<double> Psolar;
    vector<double> Pwind;
    vector<double> Pebatt;
    vector<double> Pgen;
    vector<double> Pcalc;
    double socB5 = 0.0;
    double socB7 = 0.0;
    int iterations = 0;
};

static inline double deg(double rad) {
    return rad * 180.0 / PI;
}

static inline double clampd(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

vector<vector<complex<double>>> buildYbus(int n, const vector<Branch>& branches) {
    vector<vector<complex<double>>> Y(n, vector<complex<double>>(n, complex<double>(0.0, 0.0)));
    for (const auto& br : branches) {
        complex<double> y = 1.0 / complex<double>(br.r, br.x);
        Y[br.i][br.i] += y;
        Y[br.j][br.j] += y;
        Y[br.i][br.j] -= y;
        Y[br.j][br.i] -= y;
    }
    return Y;
}

void calcPowerInjections(
    const vector<vector<complex<double>>>& Y,
    const vector<double>& Vm,
    const vector<double>& Va,
    vector<double>& P,
    vector<double>& Q
) {
    const int n = static_cast<int>(Vm.size());
    P.assign(n, 0.0);
    Q.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const double G = Y[i][k].real();
            const double B = Y[i][k].imag();
            const double a = Va[i] - Va[k];
            P[i] += Vm[i] * Vm[k] * (G * std::cos(a) + B * std::sin(a));
            Q[i] += Vm[i] * Vm[k] * (G * std::sin(a) - B * std::cos(a));
        }
    }
}

bool solveLinearSystem(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    const int n = static_cast<int>(b.size());

    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::abs(A[col][col]);
        for (int r = col + 1; r < n; ++r) {
            double v = std::abs(A[r][col]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (best < 1e-14) return false;

        if (piv != col) {
            std::swap(A[piv], A[col]);
            std::swap(b[piv], b[col]);
        }

        const double diag = A[col][col];
        for (int j = col; j < n; ++j) A[col][j] /= diag;
        b[col] /= diag;

        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            const double f = A[r][col];
            if (std::abs(f) < 1e-18) continue;
            for (int j = col; j < n; ++j) A[r][j] -= f * A[col][j];
            b[r] -= f * b[col];
        }
    }

    x = std::move(b);
    return true;
}

NRSolution newtonRaphson(
    const vector<vector<complex<double>>>& Y,
    const vector<double>& Pspec,
    const vector<double>& Qspec,
    int slack,
    const vector<int>& pvBuses,
    const vector<int>& pqBuses,
    const vector<double>& Vm0,
    const vector<double>& Va0,
    double tol = 1e-8,
    int maxIter = 40
) {
    const int n = static_cast<int>(Y.size());
    vector<double> Vm = Vm0;
    vector<double> Va = Va0;
    vector<double> P(n), Q(n);

    vector<int> pIdx;
    pIdx.reserve(n - 1);
    for (int i = 0; i < n; ++i) {
        if (i != slack) pIdx.push_back(i);
    }
    vector<int> qIdx = pqBuses;

    for (int it = 1; it <= maxIter; ++it) {
        calcPowerInjections(Y, Vm, Va, P, Q);

        vector<double> mismatch;
        mismatch.reserve(pIdx.size() + qIdx.size());

        double maxAbs = 0.0;
        for (int idx : pIdx) {
            double d = Pspec[idx] - P[idx];
            mismatch.push_back(d);
            maxAbs = std::max(maxAbs, std::abs(d));
        }
        for (int idx : qIdx) {
            double d = Qspec[idx] - Q[idx];
            mismatch.push_back(d);
            maxAbs = std::max(maxAbs, std::abs(d));
        }

        if (maxAbs < tol) {
            NRSolution out;
            out.converged = true;
            out.iterations = it;
            out.Vm = std::move(Vm);
            out.Va = std::move(Va);
            out.Pcalc = std::move(P);
            out.Qcalc = std::move(Q);
            return out;
        }

        const int npvpq = static_cast<int>(pIdx.size());
        const int npq = static_cast<int>(qIdx.size());
        const int dim = npvpq + npq;

        vector<vector<double>> J(dim, vector<double>(dim, 0.0));

        for (int r = 0; r < npvpq; ++r) {
            int i = pIdx[r];

            for (int c = 0; c < npvpq; ++c) {
                int k = pIdx[c];
                if (i == k) {
                    J[r][c] = -Q[i] - Y[i][i].imag() * Vm[i] * Vm[i];
                } else {
                    const double G = Y[i][k].real();
                    const double B = Y[i][k].imag();
                    const double a = Va[i] - Va[k];
                    J[r][c] = Vm[i] * Vm[k] * (G * std::sin(a) - B * std::cos(a));
                }
            }

            for (int c = 0; c < npq; ++c) {
                int k = qIdx[c];
                if (i == k) {
                    J[r][npvpq + c] = P[i] / Vm[i] + Y[i][i].real() * Vm[i];
                } else {
                    const double G = Y[i][k].real();
                    const double B = Y[i][k].imag();
                    const double a = Va[i] - Va[k];
                    J[r][npvpq + c] = Vm[i] * (G * std::cos(a) + B * std::sin(a));
                }
            }
        }

        for (int r = 0; r < npq; ++r) {
            int i = qIdx[r];

            for (int c = 0; c < npvpq; ++c) {
                int k = pIdx[c];
                if (i == k) {
                    J[npvpq + r][c] = P[i] - Y[i][i].real() * Vm[i] * Vm[i];
                } else {
                    const double G = Y[i][k].real();
                    const double B = Y[i][k].imag();
                    const double a = Va[i] - Va[k];
                    J[npvpq + r][c] = -Vm[i] * Vm[k] * (G * std::cos(a) + B * std::sin(a));
                }
            }

            for (int c = 0; c < npq; ++c) {
                int k = qIdx[c];
                if (i == k) {
                    J[npvpq + r][npvpq + c] = Q[i] / Vm[i] - Y[i][i].imag() * Vm[i];
                } else {
                    const double G = Y[i][k].real();
                    const double B = Y[i][k].imag();
                    const double a = Va[i] - Va[k];
                    J[npvpq + r][npvpq + c] = Vm[i] * (G * std::sin(a) - B * std::cos(a));
                }
            }
        }

        vector<double> dx;
        if (!solveLinearSystem(J, mismatch, dx)) {
            break;
        }

        for (int k = 0; k < npvpq; ++k) Va[pIdx[k]] += dx[k];
        for (int k = 0; k < npq; ++k) Vm[qIdx[k]] += dx[npvpq + k];

        for (int idx : qIdx) Vm[idx] = clampd(Vm[idx], 0.70, 1.30);
    }

    NRSolution fail;
    fail.converged = false;
    fail.iterations = maxIter;
    fail.Vm = std::move(Vm);
    fail.Va = std::move(Va);
    calcPowerInjections(Y, fail.Vm, fail.Va, fail.Pcalc, fail.Qcalc);
    return fail;
}

NRSolution solveWithQLimits(
    const vector<vector<complex<double>>>& Y,
    const vector<double>& Pspec,
    const vector<double>& Qspec,
    int slack,
    const vector<int>& pv,
    const vector<int>& pq,
    const std::unordered_map<int, std::pair<double, double>>& qLimits,
    const vector<double>& Vm0,
    const vector<double>& Va0,
    double tol,
    int maxIter
) {
    vector<int> activePv = pv;
    vector<int> activePq = pq;
    std::sort(activePq.begin(), activePq.end());
    vector<double> Qwork = Qspec;
    vector<int> qLimited;

    for (size_t t = 0; t <= pv.size(); ++t) {
        auto sol = newtonRaphson(Y, Pspec, Qwork, slack, activePv, activePq, Vm0, Va0, tol, maxIter);
        if (!sol.converged) {
            sol.activePv = activePv;
            sol.activePq = activePq;
            sol.qLimited = qLimited;
            return sol;
        }

        vector<std::pair<int, double>> violations;
        for (int b : activePv) {
            auto it = qLimits.find(b);
            double qmin = -1e9, qmax = 1e9;
            if (it != qLimits.end()) {
                qmin = it->second.first;
                qmax = it->second.second;
            }
            if (sol.Qcalc[b] < qmin) violations.push_back({b, qmin});
            else if (sol.Qcalc[b] > qmax) violations.push_back({b, qmax});
        }

        if (violations.empty()) {
            sol.activePv = activePv;
            sol.activePq = activePq;
            sol.qLimited = qLimited;
            return sol;
        }

        for (const auto& v : violations) {
            int b = v.first;
            double qfix = v.second;
            activePv.erase(std::remove(activePv.begin(), activePv.end(), b), activePv.end());
            if (std::find(activePq.begin(), activePq.end(), b) == activePq.end()) activePq.push_back(b);
            Qwork[b] = qfix;
            qLimited.push_back(b);
        }
        std::sort(activePq.begin(), activePq.end());
    }

    auto sol = newtonRaphson(Y, Pspec, Qwork, slack, activePv, activePq, Vm0, Va0, tol, maxIter);
    sol.activePv = activePv;
    sol.activePq = activePq;
    sol.qLimited = qLimited;
    return sol;
}

double solarProfile(double hour, double cloudPhase) {
    double sun = std::sin(PI * (hour - 6.0) / 12.0);
    double clearSky = std::max(0.0, sun);
    double cloud = 0.88 + 0.17 * std::sin(2.0 * PI * hour / 2.8 + cloudPhase);
    return clampd(clearSky * cloud, 0.0, 1.05);
}

double windPowerCurve(double v, double pRated, double cutIn = 3.0, double rated = 12.0, double cutOut = 25.0) {
    if (v < cutIn || v >= cutOut) return 0.0;
    if (v < rated) {
        double x = (v - cutIn) / (rated - cutIn);
        return pRated * x * x * x;
    }
    return pRated;
}

double evChargingProfile(double hour) {
    double morning = std::exp(-0.5 * std::pow((hour - 8.0) / 1.7, 2.0));
    double evening = 1.8 * std::exp(-0.5 * std::pow((hour - 19.0) / 2.2, 2.0));
    return 0.08 + 0.55 * morning + 0.95 * evening;
}

std::pair<double, double> batteryDispatch(
    double hour,
    double soc,
    double pMax,
    double eMax,
    double dtH,
    double systemStress,
    double solarStrength
) {
    const double socMin = 0.15;
    const double socMax = 0.92;

    double eveningBoost = 1.1 * std::exp(-0.5 * std::pow((hour - 19.0) / 2.6, 2.0));
    double stressTerm = 0.9 * std::max(0.0, systemStress - 1.0);
    double chargeBias = 1.2 * std::max(0.0, solarStrength - 0.45);
    double socTarget = 0.58 + 0.12 * std::sin(2.0 * PI * (hour - 3.0) / 24.0);
    double socTerm = 1.6 * (soc - socTarget);

    double pCmd = clampd(pMax * (eveningBoost + stressTerm + socTerm - chargeBias), -pMax, pMax);

    double maxDischarge = std::max(0.0, (soc - socMin) * eMax / dtH);
    double maxCharge = std::max(0.0, (socMax - soc) * eMax / dtH);
    pCmd = std::min(pCmd, maxDischarge);
    pCmd = std::max(pCmd, -maxCharge);

    double socNew = clampd(soc - (pCmd * dtH / eMax), socMin, socMax);
    return {pCmd, socNew};
}

void generateTimeSeries(int nFrames, vector<double>& hours, vector<double>& loadCurve, vector<double>& wind1, vector<double>& wind2) {
    hours.resize(nFrames);
    loadCurve.resize(nFrames);
    wind1.resize(nFrames);
    wind2.resize(nFrames);

    for (int k = 0; k < nFrames; ++k) {
        hours[k] = 24.0 * static_cast<double>(k) / static_cast<double>(nFrames);
        loadCurve[k] = clampd(
            0.86
            + 0.20 * std::sin(2.0 * PI * (hours[k] - 6.0) / 24.0)
            + 0.08 * std::sin(4.0 * PI * (hours[k] - 15.0) / 24.0),
            0.60,
            1.25
        );
    }

    wind1[0] = 8.5;
    wind2[0] = 9.8;

    uint64_t state = 88172645463325252ull;
    auto noise = [&state]() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        double u = (state % 1000000) / 1000000.0;
        return (u - 0.5) * 2.0;
    };

    for (int k = 1; k < nFrames; ++k) {
        double target1 = 8.8 + 1.8 * std::sin(2.0 * PI * hours[k] / 9.0);
        double target2 = 9.2 + 2.2 * std::sin(2.0 * PI * hours[k] / 11.0 + 0.9);
        wind1[k] = 0.88 * wind1[k - 1] + 0.12 * target1 + 0.35 * noise();
        wind2[k] = 0.86 * wind2[k - 1] + 0.14 * target2 + 0.40 * noise();
        wind1[k] = clampd(wind1[k], 2.0, 22.0);
        wind2[k] = clampd(wind2[k], 2.0, 22.0);
    }
}

bool simulateSystem(int nFrames, vector<Frame>& frames, vector<string>& busNames) {
    const int nBuses = 8;
    busNames = {
        "B1-Grid", "B2-Thermal", "B3-Thermal", "B4-UrbanLoad",
        "B5-SolarHub", "B6-WindHub", "B7-Industrial", "B8-RuralMix"
    };

    vector<Branch> branches = {
        {0, 1, 0.012, 0.085}, {0, 2, 0.016, 0.110}, {1, 3, 0.010, 0.074},
        {1, 4, 0.014, 0.091}, {2, 5, 0.012, 0.081}, {2, 6, 0.018, 0.118},
        {3, 4, 0.009, 0.063}, {4, 5, 0.013, 0.095}, {5, 6, 0.011, 0.082},
        {6, 7, 0.012, 0.089}, {3, 7, 0.017, 0.123}, {4, 7, 0.015, 0.106}
    };

    auto Y = buildYbus(nBuses, branches);
    const int slack = 0;
    vector<int> pv = {1, 2};
    vector<int> pq = {3, 4, 5, 6, 7};

    std::unordered_map<int, std::pair<double, double>> qLimits = {
        {1, {-0.50, 0.60}},
        {2, {-0.45, 0.55}},
    };

    vector<double> vmSet = {1.06, 1.03, 1.02, 1.00, 1.00, 1.00, 1.00, 1.00};
    vector<double> vmPrev = vmSet;
    vector<double> vaPrev(nBuses, 0.0);

    vector<double> baseLoadP = {0.00, 0.25, 0.22, 0.58, 0.44, 0.39, 0.52, 0.36};
    vector<double> baseLoadQ = {0.00, 0.09, 0.08, 0.22, 0.17, 0.15, 0.20, 0.13};

    const double g1Base = 1.05;
    const double g2Base = 0.85;

    const double solarB5Rated = 0.65;
    const double solarB8Rated = 0.30;
    const double windB6Rated = 0.75;
    const double windB7Rated = 0.55;

    vector<double> phase = {0.0, 0.3, 0.7, 1.4, 2.0, 2.6, 3.2, 3.9};

    std::unordered_map<int, std::pair<double, double>> batt; // bus -> (soc, emax)
    batt[4] = {0.58, 1.20};
    batt[6] = {0.52, 0.95};
    std::unordered_map<int, double> battPmax = {{4, 0.32}, {6, 0.28}};

    vector<double> hours, loadCurve, wind1, wind2;
    generateTimeSeries(nFrames, hours, loadCurve, wind1, wind2);
    const double dtH = 24.0 / static_cast<double>(nFrames);

    const double evPf = 0.96;
    const double evQFactor = std::tan(std::acos(evPf));
    std::unordered_map<int, double> evHubs = {{3, 0.30}, {6, 0.26}, {7, 0.18}};

    frames.clear();
    frames.reserve(nFrames);

    for (int k = 0; k < nFrames; ++k) {
        double hour = hours[k];

        vector<double> Pload(nBuses), Qload(nBuses);
        for (int b = 0; b < nBuses; ++b) {
            double localMod = 1.0 + 0.06 * std::sin(2.0 * PI * hour / 3.5 + phase[b]);
            double scale = clampd(loadCurve[k] * localMod, 0.55, 1.40);
            Pload[b] = baseLoadP[b] * scale;
            Qload[b] = baseLoadQ[b] * scale;
        }

        vector<double> Pev(nBuses, 0.0), Qev(nBuses, 0.0);
        double evFactor = evChargingProfile(hour);
        for (const auto& kv : evHubs) {
            int b = kv.first;
            double pmax = kv.second;
            Pev[b] = pmax * evFactor * (0.92 + 0.10 * std::sin(2.0 * PI * hour / 6.0 + b));
            Qev[b] = Pev[b] * evQFactor;
            Pload[b] += Pev[b];
            Qload[b] += Qev[b];
        }

        vector<double> Psolar(nBuses, 0.0), Pwind(nBuses, 0.0), Pbatt(nBuses, 0.0), Qbatt(nBuses, 0.0);
        Psolar[4] = solarB5Rated * solarProfile(hour, 0.2);
        Psolar[7] = solarB8Rated * solarProfile(hour, 1.1);
        Pwind[5] = windPowerCurve(wind1[k], windB6Rated);
        Pwind[6] = windPowerCurve(wind2[k], windB7Rated);

        double solarStrength = (Psolar[4] + Psolar[7]) / (solarB5Rated + solarB8Rated + 1e-9);
        double systemStress = std::accumulate(Pload.begin() + 1, Pload.end(), 0.0)
                            / (std::accumulate(baseLoadP.begin() + 1, baseLoadP.end(), 0.0) + 1e-9);

        for (auto& kv : batt) {
            int b = kv.first;
            double soc = kv.second.first;
            double emax = kv.second.second;
            auto [pCmd, socNew] = batteryDispatch(hour, soc, battPmax[b], emax, dtH, systemStress, solarStrength);
            kv.second.first = socNew;
            Pbatt[b] = pCmd;
            Qbatt[b] = 0.10 * pCmd;
        }

        vector<double> Pgen(nBuses, 0.0), Qgen(nBuses, 0.0);

        double renewTotal = 0.0;
        for (int b = 0; b < nBuses; ++b) renewTotal += Psolar[b] + Pwind[b] + Pbatt[b];

        double netDemandNoSlack = std::accumulate(Pload.begin() + 1, Pload.end(), 0.0) - renewTotal;
        double gTarget = 1.75 + 0.18 * std::sin(2.0 * PI * (hour - 1.5) / 24.0);
        double agcErr = netDemandNoSlack - gTarget;

        double g1 = clampd(g1Base + 0.45 * agcErr, 0.55, 1.45);
        double g2 = clampd(g2Base + 0.30 * agcErr, 0.40, 1.20);

        Pgen[1] = g1;
        Pgen[2] = g2;
        for (int b = 0; b < nBuses; ++b) Pgen[b] += Psolar[b] + Pwind[b] + Pbatt[b];

        Qgen[4] = 0.05 * Psolar[4];
        Qgen[7] = 0.05 * Psolar[7];
        Qgen[5] = 0.03 * Pwind[5];
        Qgen[6] = 0.03 * Pwind[6];
        for (int b = 0; b < nBuses; ++b) Qgen[b] += Qbatt[b];

        vector<double> Pspec(nBuses), Qspec(nBuses);
        for (int b = 0; b < nBuses; ++b) {
            Pspec[b] = Pgen[b] - Pload[b];
            Qspec[b] = Qgen[b] - Qload[b];
        }
        Pspec[slack] = 0.0;
        Qspec[slack] = 0.0;
        for (int b : pv) Qspec[b] = 0.0;

        vmPrev[slack] = vmSet[slack];
        vmPrev[1] = vmSet[1];
        vmPrev[2] = vmSet[2];

        auto sol = solveWithQLimits(Y, Pspec, Qspec, slack, pv, pq, qLimits, vmPrev, vaPrev, 1e-9, 45);
        if (!sol.converged) {
            std::cerr << "Frame " << k << " did not converge." << endl;
            return false;
        }

        Frame f;
        f.hour = hour;
        f.Vm = sol.Vm;
        f.VaDeg.resize(nBuses);
        for (int i = 0; i < nBuses; ++i) f.VaDeg[i] = deg(sol.Va[i]);
        f.Pload = std::move(Pload);
        f.Psolar = std::move(Psolar);
        f.Pwind = std::move(Pwind);
        f.Pebatt = std::move(Pbatt);
        f.Pgen = std::move(Pgen);
        f.Pcalc = sol.Pcalc;
        f.socB5 = batt[4].first;
        f.socB7 = batt[6].first;
        f.iterations = sol.iterations;

        frames.push_back(std::move(f));
        vmPrev = sol.Vm;
        vaPrev = sol.Va;
    }

    return true;
}

void writeCsv(const string& prefix, const vector<Frame>& frames) {
    std::ofstream summary(prefix + "_summary.csv");
    summary << "frame,hour,total_load_p,total_solar_p,total_wind_p,total_battery_p,thermal_p,slack_p,soc_b5,soc_b7,iters\n";

    std::ofstream bus(prefix + "_bus.csv");
    bus << "frame,hour,bus,vm,va_deg,p_load,p_gen,p_solar,p_wind,p_battery\n";

    for (size_t k = 0; k < frames.size(); ++k) {
        const auto& f = frames[k];
        double totalLoad = std::accumulate(f.Pload.begin(), f.Pload.end(), 0.0);
        double totalSolar = std::accumulate(f.Psolar.begin(), f.Psolar.end(), 0.0);
        double totalWind = std::accumulate(f.Pwind.begin(), f.Pwind.end(), 0.0);
        double totalBatt = std::accumulate(f.Pebatt.begin(), f.Pebatt.end(), 0.0);
        double thermal = f.Pgen[1] + f.Pgen[2];

        summary
            << k << ',' << f.hour << ',' << totalLoad << ',' << totalSolar << ','
            << totalWind << ',' << totalBatt << ',' << thermal << ',' << f.Pcalc[0] << ','
            << f.socB5 << ',' << f.socB7 << ',' << f.iterations << '\n';

        for (size_t b = 0; b < f.Vm.size(); ++b) {
            bus
                << k << ',' << f.hour << ',' << (b + 1) << ',' << f.Vm[b] << ',' << f.VaDeg[b] << ','
                << f.Pload[b] << ',' << f.Pgen[b] << ',' << f.Psolar[b] << ',' << f.Pwind[b] << ',' << f.Pebatt[b]
                << '\n';
        }
    }
}

void writeHtmlVisualization(const string& htmlPath, const vector<Frame>& frames, const vector<string>& busNames) {
    std::ofstream out(htmlPath);
    if (!out) return;

    out << "<!doctype html><html><head><meta charset='utf-8'>"
        << "<title>Load Flow Visual Studio</title>"
        << "<style>"
        << ":root{--bg:#070d1a;--panel:#0f1729;--panel2:#101a30;--line:#253451;--txt:#d9e3f0;--muted:#95a5bf;}"
        << "*{box-sizing:border-box}body{margin:0;background:radial-gradient(1200px 700px at 10% 0%,#152240,#070d1a);color:var(--txt);font-family:Segoe UI,Arial,sans-serif;}"
        << ".wrap{display:grid;grid-template-columns:2.1fr 1fr;gap:12px;padding:12px;height:100vh}"
        << ".card{background:linear-gradient(180deg,var(--panel),var(--panel2));border:1px solid #1f2d46;border-radius:12px;padding:10px;box-shadow:0 10px 28px rgba(0,0,0,.35)}"
        << ".title{font-size:14px;font-weight:700;letter-spacing:.2px;margin:0 0 8px 0;color:#e8eef9}"
        << "canvas{width:100%;display:block;border-radius:10px;background:#081224;border:1px solid #1b2a44}"
        << "#net{height:74vh}#chart{height:48vh}.statgrid{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-top:8px}"
        << ".pill{background:#0b1426;border:1px solid #223251;border-radius:8px;padding:8px}.k{font-size:11px;color:var(--muted)}.v{font-size:18px;font-weight:700}"
        << ".ctrl{display:grid;grid-template-columns:auto 1fr auto auto;gap:8px;align-items:center;margin-top:8px}"
        << "button,input{background:#111d33;color:var(--txt);border:1px solid #2a3d5d;border-radius:8px;padding:7px 10px}"
        << ".legend{font-size:12px;color:#b9c7db;margin-top:8px;line-height:1.6}.small{font-size:12px;color:#9fb0c9;margin-top:6px}"
        << "</style></head><body>";

    out << "<div class='wrap'>"
        << "<div class='card'><h3 class='title'>Dynamic 8-Bus Load Flow Animation</h3><canvas id='net' width='1180' height='760'></canvas>"
        << "<div class='ctrl'><button id='play'>Pause</button><input id='frame' type='range' min='0' max='0' value='0'>"
        << "<input id='spd' type='range' min='16' max='220' value='80'><span id='info'></span></div>"
        << "<div class='small'>Line glow, direction particles, and voltage heat are frame-by-frame from C++ simulation output.</div></div>"
        << "<div class='card'><h3 class='title'>Power and State Trends</h3><canvas id='chart' width='560' height='520'></canvas>"
        << "<div class='legend'>"
        << "<span style='color:#ff5a5f'>Load</span> | "
        << "<span style='color:#60a5fa'>Thermal</span> | "
        << "<span style='color:#f59e0b'>Solar</span> | "
        << "<span style='color:#34d399'>Wind</span> | "
        << "<span style='color:#22d3ee'>Battery</span></div>"
        << "<div class='statgrid'>"
        << "<div class='pill'><div class='k'>Total Load</div><div class='v' id='mLoad'>-</div></div>"
        << "<div class='pill'><div class='k'>Thermal</div><div class='v' id='mTherm'>-</div></div>"
        << "<div class='pill'><div class='k'>Renewables</div><div class='v' id='mRen'>-</div></div>"
        << "<div class='pill'><div class='k'>Battery SOC</div><div class='v' id='mSoc'>-</div></div>"
        << "</div></div></div>";

    out << "<script>";
    out << "const busNames=[";
    for (size_t i = 0; i < busNames.size(); ++i) {
        if (i) out << ",";
        out << "'" << busNames[i] << "'";
    }
    out << "];";

    out << "const frames=[";
    for (size_t k = 0; k < frames.size(); ++k) {
        const auto& f = frames[k];
        if (k) out << ",";
        auto vec = [&](const vector<double>& v) {
            out << "[";
            for (size_t i = 0; i < v.size(); ++i) {
                if (i) out << ",";
                out << v[i];
            }
            out << "]";
        };
        out << "{h:" << f.hour << ",vm:"; vec(f.Vm);
        out << ",va:"; vec(f.VaDeg);
        out << ",pl:"; vec(f.Pload);
        out << ",ps:"; vec(f.Psolar);
        out << ",pw:"; vec(f.Pwind);
        out << ",pb:"; vec(f.Pebatt);
        out << ",pg:"; vec(f.Pgen);
        out << ",soc5:" << f.socB5 << ",soc7:" << f.socB7 << "}";
    }
    out << "];";

    out << "const pos=[[120,370],[300,210],[300,530],[530,160],[670,290],[680,460],[510,600],[940,360]];";
    out << "const edges=[[0,1],[0,2],[1,3],[1,4],[2,5],[2,6],[3,4],[4,5],[5,6],[6,7],[3,7],[4,7]];";
    out << "const net=document.getElementById('net'), chart=document.getElementById('chart');"
        << "const nctx=net.getContext('2d'), cctx=chart.getContext('2d');"
        << "const info=document.getElementById('info'), playBtn=document.getElementById('play'), spd=document.getElementById('spd'),frameSlider=document.getElementById('frame');"
        << "const mLoad=document.getElementById('mLoad'),mTherm=document.getElementById('mTherm'),mRen=document.getElementById('mRen'),mSoc=document.getElementById('mSoc');"
        << "let i=0, run=true, timer=0; frameSlider.max=String(frames.length-1);"
        << "playBtn.onclick=()=>{run=!run;playBtn.textContent=run?'Pause':'Play';};"
        << "frameSlider.oninput=()=>{i=Number(frameSlider.value);};";

    out << "function colorVm(v){const t=Math.max(0,Math.min(1,(v-0.9)/0.2));const r=Math.floor(255*(0.2+0.8*t));const g=Math.floor(90+120*(1-Math.abs(t-0.5)*2));const b=Math.floor(255*(1-t));return `rgb(${r},${g},${b})`;}";
    out << "const totals=frames.map(f=>({"
        << "load:f.pl.reduce((a,b)=>a+b,0),"
        << "solar:f.ps.reduce((a,b)=>a+b,0),"
        << "wind:f.pw.reduce((a,b)=>a+b,0),"
        << "batt:f.pb.reduce((a,b)=>a+b,0),"
        << "thermal:f.pg[1]+f.pg[2]"
        << "}));";
    out << "const maxY=Math.max(...totals.map(t=>Math.max(t.load,t.thermal,t.solar,t.wind,Math.abs(t.batt))))*1.25;";

    out << "function flowProxy(f,a,b){const d=(f.va[a]-f.va[b]);const v=(f.vm[a]+f.vm[b])*0.5;return d*0.12*v;}"
        << "function drawArrow(x1,y1,x2,y2,color,s){const dx=x2-x1,dy=y2-y1,L=Math.hypot(dx,dy)||1,ux=dx/L,uy=dy/L;"
        << "const px=x1+ux*s,py=y1+uy*s; nctx.fillStyle=color; nctx.beginPath();"
        << "nctx.moveTo(px+ux*8,py+uy*8); nctx.lineTo(px-uy*4,py+ux*4); nctx.lineTo(px+uy*4,py-ux*4); nctx.closePath(); nctx.fill();}"
        << "function drawNet(f,t){"
        << "nctx.clearRect(0,0,net.width,net.height);"
        << "const bg=nctx.createLinearGradient(0,0,0,net.height);bg.addColorStop(0,'#0c162c');bg.addColorStop(1,'#091225');nctx.fillStyle=bg;nctx.fillRect(0,0,net.width,net.height);"
        << "for(const [a,b] of edges){const x1=pos[a][0],y1=pos[a][1],x2=pos[b][0],y2=pos[b][1];const p=flowProxy(f,a,b);"
        << "const mag=Math.min(1,Math.abs(p)/0.3);const col=p>=0?'rgba(255,102,115,':'rgba(66,153,225,';"
        << "nctx.lineWidth=2+5*mag; nctx.strokeStyle=col+(0.35+0.55*mag)+')';"
        << "nctx.beginPath();nctx.moveTo(x1,y1);nctx.lineTo(x2,y2);nctx.stroke();"
        << "nctx.strokeStyle=col+'0.22)';nctx.lineWidth=10*mag; nctx.stroke();"
        << "const L=Math.hypot(x2-x1,y2-y1)||1; const s=((t*0.12)%1)*L; if(p<0) drawArrow(x2,y2,x1,y1,col+'0.95)',s); else drawArrow(x1,y1,x2,y2,col+'0.95)',s);"
        << "}"
        << "for(let b=0;b<pos.length;b++){ const inj=f.pg[b]-f.pl[b]; const r=22+Math.min(12,Math.abs(inj)*16);"
        << "const x=pos[b][0],y=pos[b][1],c=colorVm(f.vm[b]);"
        << "nctx.shadowBlur=20; nctx.shadowColor=c; nctx.fillStyle=c; nctx.beginPath(); nctx.arc(x,y,r,0,Math.PI*2); nctx.fill();"
        << "nctx.shadowBlur=0; nctx.strokeStyle='rgba(255,255,255,0.9)'; nctx.lineWidth=1.2; nctx.stroke();"
        << "nctx.fillStyle='#06111f'; nctx.font='700 13px Segoe UI'; nctx.textAlign='center'; nctx.fillText(String(b+1),x,y+4);"
        << "nctx.fillStyle='#dbe8ff'; nctx.font='12px Segoe UI'; nctx.fillText(`${busNames[b]}  |V| ${f.vm[b].toFixed(3)}  a ${f.va[b].toFixed(2)} deg`,x,y+r+18);"
        << "}"
        << "info.textContent=`Frame ${i+1}/${frames.length}   t=${f.h.toFixed(2)} h`;"
        << "}";

    out << "function drawLine(arr,color,idx,w,h,pad){cctx.strokeStyle=color;cctx.lineWidth=2.2;cctx.beginPath();"
        << "for(let k=0;k<arr.length;k++){const x=pad+k*(w-2*pad)/(arr.length-1);const y=h-pad-(arr[k]/maxY)*(h-2*pad);if(k===0)cctx.moveTo(x,y);else cctx.lineTo(x,y);}cctx.stroke();"
        << "const x=pad+idx*(w-2*pad)/(arr.length-1);cctx.strokeStyle='rgba(226,232,240,.55)';cctx.beginPath();cctx.moveTo(x,pad);cctx.lineTo(x,h-pad);cctx.stroke();}";

    out << "function drawChart(idx){const w=chart.width,h=chart.height,pad=42;cctx.clearRect(0,0,w,h);"
        << "const g=cctx.createLinearGradient(0,0,0,h);g.addColorStop(0,'#0d1830');g.addColorStop(1,'#091222');cctx.fillStyle=g;cctx.fillRect(0,0,w,h);"
        << "cctx.strokeStyle='rgba(85,107,143,.35)';for(let gy=0;gy<=6;gy++){const y=pad+gy*(h-2*pad)/6;cctx.beginPath();cctx.moveTo(pad,y);cctx.lineTo(w-pad,y);cctx.stroke();}"
        << "drawLine(totals.map(t=>t.load),'#ff5a5f',idx,w,h,pad);"
        << "drawLine(totals.map(t=>t.thermal),'#60a5fa',idx,w,h,pad);"
        << "drawLine(totals.map(t=>t.solar),'#f59e0b',idx,w,h,pad);"
        << "drawLine(totals.map(t=>t.wind),'#34d399',idx,w,h,pad);"
        << "drawLine(totals.map(t=>Math.abs(t.batt)),'#22d3ee',idx,w,h,pad);"
        << "cctx.fillStyle='#d8e4f5';cctx.font='13px Segoe UI';cctx.fillText('Power Trends (pu)',pad,24);"
        << "}";

    out << "function updateMetrics(f){const load=f.pl.reduce((a,b)=>a+b,0),solar=f.ps.reduce((a,b)=>a+b,0),wind=f.pw.reduce((a,b)=>a+b,0),batt=f.pb.reduce((a,b)=>a+b,0),therm=f.pg[1]+f.pg[2];"
        << "mLoad.textContent=load.toFixed(2)+' pu'; mTherm.textContent=therm.toFixed(2)+' pu';"
        << "mRen.textContent=(solar+wind).toFixed(2)+' pu'; mSoc.textContent=((0.5*(f.soc5+f.soc7))*100).toFixed(1)+' %';}"
        << "function tick(ts){if(!timer)timer=ts;const step=Number(spd.value);if(run&&ts-timer>step){i=(i+1)%frames.length;frameSlider.value=String(i);timer=ts;}"
        << "const f=frames[i];drawNet(f,ts);drawChart(i);updateMetrics(f);requestAnimationFrame(tick);}requestAnimationFrame(tick);";
    out << "</script></body></html>";
}

int main(int argc, char** argv) {
    int framesN = 288;
    string csvPrefix = "load_flow_cpp";

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--frames" && i + 1 < argc) {
            framesN = std::max(24, std::stoi(argv[++i]));
        } else if (arg == "--csv-prefix" && i + 1 < argc) {
            csvPrefix = argv[++i];
        }
    }

    vector<Frame> frames;
    vector<string> busNames;

    if (!simulateSystem(framesN, frames, busNames)) {
        std::cerr << "Simulation failed." << endl;
        return 1;
    }

    writeCsv(csvPrefix, frames);
    writeHtmlVisualization(csvPrefix + "_viz.html", frames, busNames);

    const auto& last = frames.back();
    cout << "Converged frames: " << frames.size() << endl;
    cout << "Last frame t=" << std::fixed << std::setprecision(2) << last.hour << " h" << endl;
    for (size_t b = 0; b < busNames.size(); ++b) {
        cout
            << "Bus " << (b + 1) << " (" << busNames[b] << ") "
            << "|V|=" << std::setprecision(4) << last.Vm[b]
            << " angle=" << std::setprecision(3) << last.VaDeg[b] << " deg"
            << endl;
    }
    cout << "CSV written: " << csvPrefix << "_summary.csv, " << csvPrefix << "_bus.csv" << endl;
    cout << "Visualization written: " << csvPrefix << "_viz.html" << endl;

    return 0;
}
