#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include "radar/minmax_pyramid.hpp"

using namespace radar;

static RBuffer random_trace(std::size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u(-50.0f, 50.0f);
    RBuffer d(n);
    for (auto& v : d) v = u(rng);
    return d;
}

static void brute(const RBuffer& d, std::size_t i0, std::size_t i1, std::size_t columns,
                  std::vector<float>& mn, std::vector<float>& mx) {
    i1 = std::min(i1, d.size());
    const std::size_t span = i1 - i0, cols = std::min(columns, span);
    mn.assign(cols, 0.0f); mx.assign(cols, 0.0f);
    for (std::size_t c = 0; c < cols; ++c) {
        const std::size_t a = i0 + span * c / cols, b = i0 + span * (c + 1) / cols;
        float lo = d[a], hi = d[a];
        for (std::size_t i = a; i < b; ++i) { lo = std::min(lo, d[i]); hi = std::max(hi, d[i]); }
        mn[c] = lo; mx[c] = hi;
    }
}

TEST(Pyramid, MatchesBruteForceOnManyWindows) {
    const std::size_t n = 100003;
    RBuffer d = random_trace(n, 1);
    MinMaxPyramid pyr(d);
    struct W { std::size_t i0, i1, cols; };
    const W ws[] = {{0, n, 1}, {0, n, 7}, {0, n, 1000}, {0, n, 4000}, {12345, 67890, 333},
                    {50000, 50017, 5}, {5, 20, 100}, {99990, 200000, 3}, {0, 16, 16}, {1, 17, 2},
                    {3, 100000, 997}, {16, 4112, 16}};
    for (const auto& w : ws) {
        std::vector<float> mn, mx, bmn, bmx;
        const std::size_t cols = pyr.query(w.i0, w.i1, w.cols, mn, mx);
        brute(d, w.i0, w.i1, w.cols, bmn, bmx);
        ASSERT_EQ(cols, bmn.size()) << w.i0 << " " << w.i1 << " " << w.cols;
        for (std::size_t c = 0; c < cols; ++c) {
            ASSERT_EQ(mn[c], bmn[c]) << "min col " << c << " window " << w.i0 << ".." << w.i1 << " x" << w.cols;
            ASSERT_EQ(mx[c], bmx[c]) << "max col " << c << " window " << w.i0 << ".." << w.i1 << " x" << w.cols;
        }
    }
}

TEST(Pyramid, SmallAndDegenerateInputs) {
    MinMaxPyramid empty;
    std::vector<float> mn, mx;
    EXPECT_EQ(empty.query(0, 10, 5, mn, mx), 0u);

    RBuffer small = random_trace(10, 2);   // below one block: no levels, raw path
    MinMaxPyramid p(small);
    EXPECT_EQ(p.query(0, 10, 10, mn, mx), 10u);
    for (std::size_t i = 0; i < 10; ++i) { EXPECT_EQ(mn[i], small[i]); EXPECT_EQ(mx[i], small[i]); }

    EXPECT_EQ(p.query(7, 7, 4, mn, mx), 0u);
    EXPECT_EQ(p.query(8, 3, 4, mn, mx), 0u);
    EXPECT_EQ(p.query(0, 10, 0, mn, mx), 0u);
}

TEST(Pyramid, NarrowSpikesSurviveHeavyDecimation) {
    const std::size_t n = 100003;
    RBuffer d(n, 0.0f);
    d[54321] = 99.0f;
    d[70001] = -80.0f;
    MinMaxPyramid pyr(std::move(d));
    std::vector<float> mn, mx;
    ASSERT_EQ(pyr.query(0, n, 100, mn, mx), 100u);
    int up = 0, down = 0;
    for (std::size_t c = 0; c < 100; ++c) {
        if (mx[c] == 99.0f) ++up;
        if (mn[c] == -80.0f) ++down;
    }
    EXPECT_EQ(up, 1);
    EXPECT_EQ(down, 1);
}

TEST(Pyramid, ScalesToTenMillionSamples) {
    const std::size_t n = 10'000'000;
    RBuffer d = random_trace(n, 3);
    d[7'654'321] = 500.0f;
    const auto t0 = std::chrono::steady_clock::now();
    MinMaxPyramid pyr(std::move(d));
    const auto t1 = std::chrono::steady_clock::now();
    std::vector<float> mn, mx;
    float best = -1e9f;
    for (int rep = 0; rep < 100; ++rep) {
        pyr.query(static_cast<std::size_t>(rep) * 1000, n - static_cast<std::size_t>(rep) * 1000, 2000, mn, mx);
        for (float v : mx) best = std::max(best, v);
    }
    const auto t2 = std::chrono::steady_clock::now();
    const double build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double query_ms = std::chrono::duration<double, std::milli>(t2 - t1).count() / 100.0;
    std::printf("[ INFO ] build %.1f ms, query %.3f ms per 2000-column frame\n", build_ms, query_ms);
    EXPECT_EQ(best, 500.0f);
    EXPECT_LT(query_ms, 5.0);
}

TEST(Pyramid, BorrowedStridedViewMatchesBruteForceAndSurvivesMove) {
    const std::size_t n = 50011;
    RBuffer inter = random_trace(2 * n, 9);
    RBuffer iv(n), qv(n);
    for (std::size_t i = 0; i < n; ++i) { iv[i] = inter[2 * i]; qv[i] = inter[2 * i + 1]; }
    MinMaxPyramid pi(inter.data(), n, 2), pq(inter.data() + 1, n, 2);
    MinMaxPyramid pi_moved(std::move(pi));
    struct W { std::size_t i0, i1, cols; };
    const W ws[] = {{0, n, 500}, {1234, 40000, 77}, {0, n, 1}, {100, 140, 40}};
    for (const auto& w : ws) {
        std::vector<float> mn, mx, bmn, bmx;
        ASSERT_EQ(pi_moved.query(w.i0, w.i1, w.cols, mn, mx), std::min(w.cols, w.i1 - w.i0));
        brute(iv, w.i0, w.i1, w.cols, bmn, bmx);
        for (std::size_t c = 0; c < mn.size(); ++c) { ASSERT_EQ(mn[c], bmn[c]); ASSERT_EQ(mx[c], bmx[c]); }
        pq.query(w.i0, w.i1, w.cols, mn, mx);
        brute(qv, w.i0, w.i1, w.cols, bmn, bmx);
        for (std::size_t c = 0; c < mn.size(); ++c) { ASSERT_EQ(mn[c], bmn[c]); ASSERT_EQ(mx[c], bmx[c]); }
    }
}

TEST(Pyramid, OwningPyramidSurvivesMove) {
    RBuffer d = random_trace(1000, 4);
    MinMaxPyramid a(d);
    MinMaxPyramid b(std::move(a));
    std::vector<float> mn, mx, bmn, bmx;
    b.query(0, 1000, 50, mn, mx);
    brute(d, 0, 1000, 50, bmn, bmx);
    for (std::size_t c = 0; c < mn.size(); ++c) { EXPECT_EQ(mn[c], bmn[c]); EXPECT_EQ(mx[c], bmx[c]); }
}
