#pragma once

#include <string>
#include <type_traits>

#include "VCL2/vectorclass.h"

/*
 * Filter coefficients coef_lf and coef_hf taken from BBC PH-2071 (Weston 3 Field Deinterlacer).
 * Used when there is spatial and temporal interpolation.
 * Filter coefficients coef_sp are used when there is spatial interpolation only.
 * Adjusted for matching visual sharpness impression of spatial and temporal interpolation.
 */
static constexpr uint16_t coef_lf[2] = { 4309, 213 };
static constexpr uint16_t coef_hf[3] = { 5570, 3801, 1016 };
static constexpr uint16_t coef_sp[2] = { 5077, 981 };
static constexpr float coef_lf_f[2] = { 4309 / 8192.0f, 213 / 8192.0f };
static constexpr float coef_hf_f[3] = { 5570 / 8192.0f, 3801 / 8192.0f, 1016 / 8192.0f };
static constexpr float coef_sp_f[2] = { 5077 / 8192.0f, 981 / 8192.0f };

template<typename pixel_t, bool spat>
void filterEdge_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template<typename pixel_t>
void filterLine_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;

template<typename pixel_t, bool spat>
void filterEdge_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template<typename pixel_t>
void filterLine_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;

template<typename pixel_t, bool spat>
void filterEdge_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template<typename pixel_t>
void filterLine_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
