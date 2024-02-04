#include <algorithm>
#include <array>
#include <iostream>

#include "catch2/catch_test_macros.hpp"
#include "mini/mini.h"

namespace mini {
namespace test {

constexpr auto kTestSeq =
    "AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC"
    "AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC";

constexpr auto kTestSeqLen = 140U;

constexpr auto kKMerLen = 29U;
constexpr auto kKMerLenLong = 58U;
constexpr auto kWinLen = 9U;

constexpr auto kTestMinimizers = std::array<KMer, 23U>{
    KMer(646286436247489772, 4, 0),  KMer(4800021203147078872, 7, 0),
    KMer(3489739751552373074, 15, 0), KMer(181256282556872170, 21, 0),
    KMer(301194077553889042, 28, 0), KMer(1984613911990073889, 29, 0),
    KMer(3969227673656292348, 30, 0), KMer(5164955381998653605, 39, 0),
    KMer(3604134460204853669, 41, 0), KMer(3579537055126917398, 46, 0),
    KMer(2297258160501326761, 49, 1), KMer(3003143238680945871, 57, 0),
    KMer(1972518964064674621, 63, 1), KMer(1348383673427691446, 68, 0),
    KMer(646286436247489772, 74, 0),  KMer(4800021203147078872, 77, 0),
    KMer(3489739751552373074, 85, 0), KMer(181256282556872170, 91, 0),
    KMer(301194077553889042, 98, 0),  KMer(1984613911990073889, 99, 0),
    KMer(3969227673656292348, 100, 0), KMer(5164955381998653605, 109, 0),
    KMer(3604134460204853669, 111, 0)};

constexpr auto kTestMinimizersLong = std::array<KMer, 14U>{
    KMer(916718422212990171, 6, 0),  KMer(4626220760283687045, 11, 0),
    KMer(2591599766356106169, 17, 0), KMer(176679790298658678, 24, 0),
    KMer(185623040052900602, 27, 0), KMer(268886174752216349, 35, 0),
    KMer(2750856359964708494, 38, 0), KMer(1152133110661816871, 47, 0),
    KMer(1255617870259713102, 53, 0), KMer(749452792572590264, 57, 0),
    KMer(144082624022435439, 60, 0), KMer(1325561104408984909, 62, 0),
    KMer(2434044891487640860, 69, 0), KMer(916718422212990171, 76, 0)};

}  // namespace test
}  // namespace mini

TEST_CASE("Minimize", "[minimize]") {
  auto const minimizers = mini::Minimize(
      mini::test::kTestSeq, mini::test::kKMerLen, mini::test::kWinLen);

  REQUIRE(minimizers.size() == mini::test::kTestMinimizers.size());
  for (auto i = 0U; i < minimizers.size(); ++i) {
    CHECK(minimizers[i].value() == mini::test::kTestMinimizers[i].value());
    CHECK(minimizers[i].position() ==
          mini::test::kTestMinimizers[i].position());
    CHECK(minimizers[i].strand() == mini::test::kTestMinimizers[i].strand());
  }
}

TEST_CASE("MinimizeLong", "[minimize]") {
  auto const minimizers = mini::Minimize(
      mini::test::kTestSeq, mini::test::kKMerLenLong, mini::test::kWinLen);

  REQUIRE(minimizers.size() == mini::test::kTestMinimizersLong.size());
  for (auto i = 0U; i < minimizers.size(); ++i) {
    CHECK(minimizers[i].value() == mini::test::kTestMinimizersLong[i].value());
    CHECK(minimizers[i].position() ==
          mini::test::kTestMinimizersLong[i].position());
    CHECK(minimizers[i].strand() == mini::test::kTestMinimizersLong[i].strand());
  }
}
